function obj = doAutophase( obj)
%function obj = doAutophase( obj)
% Stoch and Olejniczak based correction
%   view by view missing point and phase estimation
%   see pmid_15705522
% CAC 201210
% Copyright Champaign Imaging LLC
% Patent Pending
% This work is licensed under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
% https://creativecommons.org/licenses/by-nc-nd/4.0/

%DEBUG_FLAG = obj.FLAGS.DEBUG;
DEBUG_FLAG = 3;

if ( DEBUG_FLAG >=3 )
    fprintf( '\n%s DEBUG ====================\n', mfilename( 'fullpath'));
end

% reply = input( 'Do you want to autophase and estimate missing points? Y/N [Y]:','s');
% if isempty( reply)
%     reply = 'Y';
% end
% if reply == 'Y'
%     % Phase all fids
% else
%     return;
% end

if ( DEBUG_FLAG >= 2 )
    tic;    fprintf( 'Processing projection info for Stoch/Olejniczak autophasing, ');
end

% temporary, CAC 210502 ***
obj.switches.est_extra = 0;     % 1 for actual data, 0 for simulations?
obj.switches.correct_mp = 1;
obj.switches.apply_phi = 0;

% Read Switches and variables
est_extra = obj.switches.est_extra;     % extra points to estimate
correct_mp = obj.switches.correct_mp;   % DEFAULT do the correction
apply_phi = obj.switches.apply_phi;     % apply the phase correction to the fids if true

startPath = obj.startPath;
savePath = obj.savePath;
%N_views = obj.views.Nviews;
%N_views_nav = obj.pars.N_views_nav;
%N_views_s = obj.pars.N_views_s;
rcvr_gate = obj.pars.rcvr_gate;

%Local switches and variables
idx_display = 1;    % view for display graphs, not working for all plots??? *** CAC 210406
ch_display = 1;     % channel for display graphs
fixed_scale = 0;    % fixed plot scales for simulation or normalized data
center_proj = 1;    % center projections between baselines at ends for phase and mp calculation

% sizes and fix for channels = 1, volumes = 1;
sz_fid = size( obj.fid);
ndims_fid = ndims( obj.fid);
if ndims( obj.fid) < 3
    sz_fid(3) = 1;
end
if ndims_fid < 4
    sz_fid(4) = 1;
end
np = sz_fid(1); nv = sz_fid(2); ns = sz_fid(3); nc = sz_fid(4);
% flatten to views
total_views = nv*ns*nc;


% Stoch and Olejniczak based correction
%   see pmid_15705522
%   CAC 200910

N = np;     % total points in fid
if rcvr_gate <= 0
    Nd = abs( rcvr_gate) + est_extra;	% missing ponts at beginning of fid
    %Nd = 3; % for testing and consistancy? *** CAC 201130
    Nb = N*7/16;    % number of baseline points at edges of total FOV, Nb/2 on each side
elseif rcvr_gate > 0
    Nd = 3; % for old datasets? *** CAC 201130
    Nb = N*5/16;    % number of baseline points at edges of total FOV, Nb/2 on each side
end

Nd_text = sprintf( '        Nd  =  %g', Nd);
Nb_text = sprintf( '        Nb  =  %g', Nb);

% put debug information into obj
if ( DEBUG_FLAG >= 3 )
    obj.debug.doAutophase.N = N;
    obj.debug.doAutophase.Nd = Nd;
    obj.debug.doAutophase.Nb = Nb;
    %obj.debug.doAutophase.Nc = Nc;
end

if ( DEBUG_FLAG >=3 )
    fprintf( '\n%s DEBUG %d missing fid points out of %d, with %d baseline points.\n', mfilename, Nd, N, Nb);
end

% Generate the A matrix Stoch eq. (12)
Nc = 1 + (N/2); % this could be put in the parfor loop, per view...., *** CAC 201210
%p_dist = fftshift( fft( fid_parfor(:, view_idx)));
%Nc = ( abs( p_dist)' * (1:N)' )  / sum( abs( p_dist) );
A = zeros( 2*Nd, 2*Nd); % this does not change for eq. (15)
alpha = zeros( Nb, 2*Nd);  alpha_star = alpha;
const = 2*pi/N;
nu_idxs = (1:N) - N/2 -1; zerofreqidx = 1 + N/2; zerofreqval = nu_idxs( zerofreqidx); % debug centering
nu_idxs_shift = Nb/2 - round( Nc-(N/2)) + 1;
nu_idxs = circshift( nu_idxs, + nu_idxs_shift);
nu_idxs = nu_idxs(1:Nb);
nu_idxs = nu_idxs';
for m_idx = 1:(2*Nd)
    m_idx_round = floor( (m_idx-1)/2);
    coson_switch_m = mod( m_idx, 2); sinon_switch_m =  (1 - coson_switch_m);
    coson_switch_n = mod( m_idx, 2); sinon_switch_n =  (1 - coson_switch_n); % why are these the same sign? It certainly works better than opposite sign *** CAC 201212
    alpha_star(:, m_idx) = coson_switch_m*cos( const*(m_idx_round) .* nu_idxs) + sinon_switch_m*sin( const*(m_idx_round) .* nu_idxs);
    alpha(:, m_idx)      = coson_switch_n*cos( const*(m_idx_round) .* nu_idxs) + sinon_switch_n*sin( const*(m_idx_round) .* nu_idxs);
end
if Nd > 0
    A  = alpha_star' * alpha;
    A(2, 2) = 1; % can't estimate the dc imaginary component with A alone due to symmetry, see Stoch (A.16)
    A_cond = cond( A);
else
    A = 1;
end
A_cond = cond( A);
A_cond_text = sprintf( '        A matrix cond  =  %g', A_cond);
cond_color = 'black';
if A_cond > 1000;
    fprintf( '\n'); warning( '"A" matrix conditon number is %g',  A_cond);
    cond_color = '#D95319'; % orage for warning
end

if ( DEBUG_FLAG >=3 )
    fprintf( '\n%s DEBUG A matrix condition number is  %g\n', mfilename, A_cond);
end

% Do SVD pseudoinverse
pinv_tol = 1.0e-5;
A_pinv = pinv( A, pinv_tol);
A_pinv_cond = cond( A_pinv);

% put debug information into obj
if ( DEBUG_FLAG >= 3 )
    obj.debug.doAutophase.alpha = alpha;
    obj.debug.doAutophase.alpha_star = alpha_star;
    obj.debug.doAutophase.A = A;
    obj.debug.doAutophase.A_cond = A_cond;
    obj.debug.doAutophase.A_pinv = A_pinv;
    obj.debug.doAutophase.pinv_tol = pinv_tol;
    obj.debug.doAutophase.A_pinv_cond = A_pinv_cond;
end

% baseline bands (f0 at element 1) Kuethe eq. [6]0
%   prep for F vectors in parfor loop
t_vec = zeros( N, 1); t_vec(1:Nb) = 1;
t_vec_shift = Nb/2 - round( Nc-(N/2)) + 1;
t_vec = circshift( t_vec, -t_vec_shift);
%t_vec = t_vec;
T_mat = diag( t_vec);
%p_sz = size( p_dist);
%T_sz = size( T_mat);

% Generate the G matrix
G = complex( eye( Nb, Nb));
for m_idx = 1:(2*Nd)
    m_idx_round = floor( (m_idx-1)/2);
    for n_idx = 1:(2*Nd)
        n_idx_round = floor( (n_idx-1)/2);
        %mn_idxs_test( m_idx, n_idx) = complex( m_idx_round, n_idx_round); % *** debug CAC 200917
        G = G - A_pinv(m_idx, n_idx) * alpha_star(:, n_idx) * alpha_star(:, m_idx)';
    end
end

G_cond = cond( G);
        % check other quadrants? *** CAC 200919
        %phi = .25*pi; % *** debug
G_cond_text = sprintf( '        G matrix cond  =  %g', G_cond);
cond_color_G = 'black';
if G_cond > 1.0e16;
    fprintf( '\n'); warning( '"G" matrix conditon number is %g',  G_cond);
    cond_color_G = '#D95319'; % orage for warning
end

if ( DEBUG_FLAG >=3 )
    fprintf( '\n%s DEBUG G matrix condition number is  %g\n', mfilename, G_cond);
end

% put debug information into obj
if ( DEBUG_FLAG >= 3 )
    obj.debug.doAutophase.G = G;
    obj.debug.doAutophase.G_cond = G_cond;
end

if ( DEBUG_FLAG >= 2 )
    fprintf( 'parfor loop running...');
end

%% use a big flat parfor loop
fid_parfor = reshape( obj.fid, np, total_views);
phi_parfor = zeros( total_views, 1);
p_dist_shift = zeros( total_views, 1);
%fid_window = fftshift( gausswin(2*N, 6)); fid_window = fid_window(1:N);
fid_before = fid_parfor(:, idx_display);

if correct_mp == 1;
    parfor view_idx = 1:total_views;
        % zero out fid missing points and window
        fid_dist = fid_parfor(:, view_idx);
        fid_dist(1:Nd) = 0.0;
        %fid_dist = fid_dist .* fid_window; % does not seem to help, *** CAC 201212
        
        % Generate and solve p*tau^2+q*tau-p=0, for eq. (15) of Stoch
        p_dist = fftshift( fft( fid_dist));
        
        if center_proj
            % find magnitude centroid and shift, WORKING HERE *** CAC 201410
            cent1 = (1:np) * abs( p_dist)/sum( abs( p_dist)); % conservative
            %cent = (1:np) * real( p_dist)/sum( real( p_dist)); % debug
            [~, cent2] = max( abs( p_dist)); % erratic
            cent = cent1;
            p_dist_shift(view_idx) = round( cent-Nc);
            %p_dist_shift(view_idx) = 0; % debug *** CAC 210408
            p_dist = circshift( p_dist, -p_dist_shift(view_idx));
            fid_dist = ifft( ifftshift( p_dist)); % do phase ramp multiply instead? *** CAC 210409
        else
            p_dist_shift(view_idx) = 0;
        end
        
        fid_pinv = fid_dist;
        p_base = T_mat * p_dist;
        
        % Shift and clip the baseline data
        W = circshift( p_base, t_vec_shift);
        W = W(1:Nb);
        W_re = real( W);
        W_im = imag( W);
        
        % compute p and q
        p = -(G * W_re)' * W_im;
        q = (G * W_re)' * W_re - (G * W_im)' * W_im;
        
        % solve for tau
        tau_all = roots( [p q -p]);
        if p < 0 % still not always picking the correct one... *** CAC 201015
            tau = tau_all(2);
        else
            tau = tau_all(1);
        end
        
        % compute phi (independent of shift, *** CAC 210409)
        phi = atan( tau); % only working for -pi/2 phi < pi/2 ***
        phi_parfor( view_idx) = phi; % store phi for debugging, *** CAC 201211
        phi_deg = phi * (180.0/pi);
        
        % Estimate missing data, see Stoch A.14
        F_re = real( exp( i*phi) * W);
        f_re = alpha_star' * F_re;
        dm_real = -A_pinv * f_re; % check the sign
        
        % Convert to complex
        dm = complex( zeros( Nd, 1));
        for idx = 1:Nd
            re_idx = 2*idx-1;
            im_idx = 2*idx;
            dm( idx) = dm_real(re_idx) + i*dm_real(im_idx);
        end
        
        % Set up the estimated fid
        %   fid_pinv defined at top of parfor loop
        fid_pinv = exp( i*phi) * fid_pinv;
        fid_pinv(1:Nd) = dm;
        %dm = exp( i*phi) * dm;

        f_ramp = ((1:N)' - 1)/(N);
        shiftramp = exp( 2.0*i*pi*p_dist_shift(view_idx)*f_ramp);
        if apply_phi
            if real( fid_pinv(1)) < 0
                fid_pinv = -fid_pinv;
            end
        else
            shiftramp = exp( -i*phi) * shiftramp; % take out the phi correction
        end
        fid_parfor(:, view_idx) = shiftramp .* fid_pinv;
    end
else
    fprintf( '\n'); warning( 'Missing point correction is disabled, using zeros!');
    fid_parfor(1:Nd, :) = complex( single( 0.0)); % just zero out the missing points, do not correct
end

obj.fid = reshape( fid_parfor, np, nv, ns, nc);
phi_parfor = reshape( phi_parfor, nv, ns, nc);

phi_avg = sum( phi_parfor(:, :, ch_display), 'all')/numel( phi_parfor(:, :, ch_display));
phi_avg_text = sprintf( '        phi avg  =  %g deg', phi_avg*180.0/pi);

if ( DEBUG_FLAG >=3 )
    fprintf( '\n%s DEBUG phi avg  =  %g deg\n', mfilename, phi_avg*180.0/pi);
end

if ( DEBUG_FLAG >= 2 )
    fprintf( 'parfor loop done, ');
end
% put debug information into obj
if ( DEBUG_FLAG >= 3 )
    obj.debug.doAutophase.phi = phi_parfor;
    obj.debug.doAutophase.phi_avg = phi_avg;
end

if ( DEBUG_FLAG >=3 )
    figure( 'name', 'Stoch/Olejniczak Phase and Missing Points', 'Position', [10 10 1376 394]);
    subplot( 1, 2, 1);
    hold on; title( 'Phase Estimate'); axis tight; xlabel( 'View [index]'); ylabel( 'Phase [deg]');
    Ndisp_seg = min( 1024, 4*nv);
    plot(  phi_parfor(1:Ndisp_seg)*180.0/pi, 'k. ', 'LineWidth', 1);
    ylim_min = ylim; ylim_max = ylim_min(2); ylim_min = ylim_min(1);
    range = ylim_max - ylim_min; sgn_ylim = sign( ylim_max);
    if sgn_ylim > 0
        start = ylim_min;
    else
        start = ylim_max; range = -range;
    end
    text( 32, start + range*0.50, phi_avg_text)
    text( 32, start + range*0.55, Nd_text);
    text( 32, start + range*0.60, Nb_text);
    %text( 32, start + range*0.65, G_cond_text, 'color', cond_color_G);
    text( 32, start + range*0.70, A_cond_text, 'color', cond_color);
    %cd( savePath); saveas( gcf, 'phase_estimate', 'png'); cd( startPath);
    
    if Nd > 0
        %figure( 'name', 'Stoch/Olejniczak missing points', 'Position', [10 10 437 330]);
        subplot( 1, 2, 2); hold on;
        title( 'Missing Points Real and Imaginary'); axis tight;
        xlabel( 'View [index]'); ylabel( 'Signal [Arb]');
        colors = lines;
        for idx = 1:Nd
            plot( real( fid_parfor(idx, 1:Ndisp_seg)), '-', 'color', colors( idx, :),'LineWidth', 1);
            plot( imag( fid_parfor(idx, 1:Ndisp_seg)), ':', 'color', colors( idx, :),'LineWidth', 1);
        end
        legend; hold off;
        cd( savePath); saveas( gcf, 'phase_and_missing_points', 'png'); cd( startPath);
        cd( savePath); saveas( gcf, 'phase_and_missing_points', 'fig'); cd( startPath); % for ismrm simulation, *** CAC 210412
    end
end

% graphs for debugging and publication
%   "fid_before" defined above
Ndisp = 16;
if ( DEBUG_FLAG >=3 )
    %figure( 'name', 'Stoch/Olejniczak missing points/phase estimate', 'Position', [10 10 1376 778]);
    figure( 'name', 'Stoch/Olejniczak FIDs and Projections', 'Position', [10 10 1376 778]);
    subplot( 2, 2, 1);
    hold on; title( 'Raw FID'); xlabel( 'Time [index]'); ylabel( 'Signal [arb]'); axis tight;
    plot(  real( fid_before(1:Ndisp)), 'g+:', 'LineWidth', 2);
    plot(  imag( fid_before(1:Ndisp)), 'rx:', 'LineWidth', 2);
    plot(  real( fid_before(1:Nd)   ), 'b+:', 'LineWidth', 2);
    plot(  imag( fid_before(1:Nd)   ), 'mx:', 'LineWidth', 2);
    xticks( 2:2:Ndisp);
    if fixed_scale; axis( [-inf inf -.13 0.3]); end; % for ismrm simulation, *** CAC 210412
    legend( 'Re data', 'Im data', 'Re missing', 'Im missing'); grid on; %hold off;
end

p_before = fftshift( fft( fid_before));
b_left_idxs = 1:(Nb/2); b_right_idxs = (N-Nb/2+1):N;

if ( DEBUG_FLAG >=3 )
    %figure( 'name', 'Stoch/Olejniczak missing points/phase estimate', 'Position', [10 10 437 330]);
    subplot( 2, 2, 3);
    hold on; title( 'Raw Projection'); xlabel( 'Frequency [index]'); ylabel( 'Signal [arb]'); axis tight;
    plot(  abs( p_before), 'k -.', 'LineWidth', 2);
    plot( real( p_before), 'g -.', 'LineWidth', 2);
    plot( imag( p_before), 'r -.', 'LineWidth', 2);
    plot( real( p_before(b_left_idxs)), 'b -', 'LineWidth', 2);
    plot( imag( p_before(b_left_idxs)), 'm -', 'LineWidth', 2);
    plot( b_right_idxs, real( p_before(b_right_idxs)), 'b -', 'LineWidth', 2);
    plot( b_right_idxs, imag( p_before(b_right_idxs)), 'm -', 'LineWidth', 2);
    xticks( (N/8):(N/8):N);
    if fixed_scale; axis( [-inf inf -1.1 1.1]); end; % for ismrm simulation, *** CAC 210412
    %legend( 'Re object', 'Im object', 'Re baseline', 'Im baseline');grid on; %hold off;
    legend( 'Mag obj', 'Re obj', 'Im obj', 'Re baseline', 'Im baseline'); grid on; %hold off;
end

phi_est = phi_parfor(idx_display);
fid_after = exp( i*phi_est)*fid_parfor( :, idx_display);
phi_text = sprintf( '        estimated phi  =  %g deg', phi_est*180.0/pi);

if ( DEBUG_FLAG >=3 )
    %figure( 'name', 'Stoch/Olejniczak missing points/phase estimate', 'Position', [10 10 437 330]);
    subplot( 2, 2, 2);
    hold on; title( 'Estimated FID'); xlabel( 'Time [index]'); ylabel( 'Signal [arb]'); axis tight;
    plot(  real( fid_after(1:Ndisp)), 'g+:', 'LineWidth', 2);
    plot(  imag( fid_after(1:Ndisp)), 'rx:', 'LineWidth', 2);
    plot(  real( fid_after(1:Nd)   ), 'b+:', 'LineWidth', 2);
    plot(  imag( fid_after(1:Nd)   ), 'mx:', 'LineWidth', 2);
    legend( 'Re data', 'Im data', 'Re est', 'Im est'); grid on; %hold off;
    if fixed_scale; axis( [-inf inf -.13 0.3]); end; % for ismrm simulation, *** CAC 210412
    ylim_min = ylim; ylim_max = ylim_min(2); ylim_min = ylim_min(1);
    range = ylim_max - ylim_min; sgn_ylim = sign( ylim_max);
    if sgn_ylim > 0
        start = ylim_min;
    else
        start = ylim_max; range = -range;
    end
    text( Ndisp/4, start + range*0.50, phi_avg_text)
    text( Ndisp/4, start + range*0.55, Nd_text);
    text( Ndisp/4, start + range*0.60, Nb_text);
    %text( Ndisp/4, start + range*0.65, G_cond_text, 'color', cond_color_G);
    text( Ndisp/4, start + range*0.70, A_cond_text, 'color', cond_color);
    grid on; %hold off;
end

p_after = fftshift( fft( fid_after));
p_after = circshift( p_after, -p_dist_shift(idx_display));%     plot(  imag( fid_before(1:Ndisp)), 'rx:', 'LineWidth', 2);
%     plot(  real( fid_before(1:Nd)   ), 'b+:', 'LineWidth', 2);
%     plot(  imag( fid_before(1:Nd)   ), 'mx:', 'LineWidth', 2);
%     xticks( 2:2:Ndisp);

if ( DEBUG_FLAG >=3 )
    %figure( 'name', 'Stoch/Olejniczak missing points/phase estimate', 'Position', [10 10 437 330]);
    subplot( 2, 2, 4);
    hold on; title( 'Estimated Projection (centered)'); xlabel( 'Frequency [index]'); ylabel( 'Signal [arb]'); axis tight;
    plot(  abs( p_after), 'k -.', 'LineWidth', 2);
    plot( real( p_after), 'g -.', 'LineWidth', 2);
    plot( imag( p_after), 'r -.', 'LineWidth', 2);
    plot( real( p_after(b_left_idxs)), 'b -', 'LineWidth', 2);
    plot( imag( p_after(b_left_idxs)), 'm -', 'LineWidth', 2);
    plot( b_right_idxs, real( p_after(b_right_idxs)), 'b -', 'LineWidth', 2);
    plot( b_right_idxs, imag( p_after(b_right_idxs)), 'm -', 'LineWidth', 2);
    xticks( (N/8):(N/8):N);
    if fixed_scale; axis( [-inf inf -1.1 1.1]); end; % for ismrm simulation, *** CAC 210412
    legend( 'Mag obj', 'Re obj', 'Im obj', 'Re baseline', 'Im baseline'); grid on; %hold off;
    grid on; hold off;
    cd( savePath); saveas( gcf, 'fid_and_baseline', 'png'); cd( startPath);
    cd( savePath); saveas( gcf, 'fid_and_baseline', 'fig'); cd( startPath); % for ismrm simulation, *** CAC 210412
end

if ( DEBUG_FLAG >=4 )
    figure( 'name', 'p_dist_shift');
    hold on; title( 'Estimated Centering'); xlabel( 'View [index]'); ylabel( 'p-dist-shift'); axis tight;
    plot( p_dist_shift(1:1024)); % debugging cent calculation, CAC 210329 ***
end

if ( DEBUG_FLAG >=4 ) & Nd > 0
    alpha_plot = padarray( alpha, N-Nb, 'pre'); % prepare alpha
    alpha_plot = circshift( alpha_plot, Nb/2);
    %figure( 'name', 'Stoch/Olejniczak missing points/phase estimate', 'Position', [10 10 1376 778]);
    figure( 'name', 'Stoch/Olejniczak Pseudinverse Related', 'Position', [10 10 1376 778]);
    subplot( 2, 2, 1);
    hold on; title( 'alpha'); xlabel( 'Frequency [index]'); ylabel( 'Signal [arb]'); axis tight;
    plot( alpha_plot(b_left_idxs, 1), 'b -', 'LineWidth', 2);
    plot( alpha_plot(b_left_idxs, 2), 'm -', 'LineWidth', 2);
    plot( b_right_idxs, alpha_plot(b_right_idxs, 1), 'b -', 'LineWidth', 2);
    plot( b_right_idxs, alpha_plot(b_right_idxs, 2), 'm -', 'LineWidth', 2);
    if Nd > 2
        plot( alpha_plot(b_left_idxs, 3), 'b --', 'LineWidth', 2);
        plot( alpha_plot(b_left_idxs, 4), 'm --', 'LineWidth', 2);
        plot( b_right_idxs, alpha_plot(b_right_idxs, 3), 'b --', 'LineWidth', 2);
        plot( b_right_idxs, alpha_plot(b_right_idxs, 4), 'm --', 'LineWidth', 2);
    end
    if Nd> 3
        plot( alpha_plot(b_left_idxs, 5), 'b :', 'LineWidth', 2);
        plot( alpha_plot(b_left_idxs, 6), 'm :', 'LineWidth', 2);
        plot( b_right_idxs, alpha_plot(b_right_idxs, 5), 'b :', 'LineWidth', 2);
        plot( b_right_idxs, alpha_plot(b_right_idxs, 6), 'm :', 'LineWidth', 2);
    end
    if fixed_scale; axis( [-inf inf -1.1 1.1]); end; % for ismrm simulation, *** CAC 210416
    legend( 'Re', 'Im');
    grid on; hold off;
    
%     alpha_plot = padarray( alpha_star, N-Nb, 'pre');
%     alpha_plot = circshift( alpha_plot, Nb/2);
%     %figure( 'name', 'Stoch/Olejniczak missing points/phase estimate', 'Position', [10 10 1376 778]);
%     %figure( 'name', 'Stoch/Olejniczak Pseudinverse Related', 'Position', [10 10 1376 778]);
%     subplot( 2, 2, 2);
%     hold on; title( 'alpha star'); xlabel( 'Frequency [index]'); ylabel( 'Signal [arb]'); axis tight;
%     plot( alpha_plot(b_left_idxs, 1), 'b -', 'LineWidth', 2);
%     plot( alpha_plot(b_left_idxs, 2), 'm -', 'LineWidth', 2);
%     plot( b_right_idxs, alpha_plot(b_right_idxs, 1), 'b -', 'LineWidth', 2);
%     plot( b_right_idxs, alpha_plot(b_right_idxs, 2), 'm -', 'LineWidth', 2);
%     if Nd > 2
%         plot( alpha_plot(b_left_idxs, 3), 'b --', 'LineWidth', 2);
%         plot( alpha_plot(b_left_idxs, 4), 'm --', 'LineWidth', 2);
%         plot( b_right_idxs, alpha_plot(b_right_idxs, 3), 'b --', 'LineWidth', 2);
%         plot( b_right_idxs, alpha_plot(b_right_idxs, 4), 'm --', 'LineWidth', 2);
%     end
%     if Nd > 3
%         plot( alpha_plot(b_left_idxs, 5), 'b :', 'LineWidth', 2);
%         plot( alpha_plot(b_left_idxs, 6), 'm :', 'LineWidth', 2);
%         plot( b_right_idxs, alpha_plot(b_right_idxs, 5), 'b :', 'LineWidth', 2);
%         plot( b_right_idxs, alpha_plot(b_right_idxs, 6), 'm :', 'LineWidth', 2);
%     end
%     if fixed_scale; axis( [-inf inf -1.1 1.1]); end; % for ismrm simulation, *** CAC 210416
%     legend( 'Re', 'Im');
%     grid on; hold off;
end

if ( DEBUG_FLAG >= 2 )
    toc
end

return

