function obj = doKdata( obj)
% ngfnRecon.doKdata
% Grids radial k-space data to cartesian k-space
%   C. Corum 5/3/2015
%
%   CAC 190226  separate DCF calculation, gridding, image manipulation and saving.
% Copyright Champaign Imaging LLC
% This work is licensed under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
% https://creativecommons.org/licenses/by-nc-nd/4.0/

%% Initial version and debug info
DEBUG_FLAG = obj.FLAGS.DEBUG;
%DEBUG_FLAG = 4;

if ( DEBUG_FLAG >=3 )
    fprintf( '\n%s DEBUG ====================\n', mfilename( 'fullpath'));
end

% version information
version = obj.version;

SIZE_FLAG = obj.FLAGS.SIZE;

if ( DEBUG_FLAG >= 2 ); kdata_timing = tic; end     % start timer

%% Switches

%debug
%switches = obj.switches
density_type = obj.switches.density_type;

%% Set up variables and data

% paths
startPath = obj.startPath;
savePath = obj.savePath;

% parameters from healpix based vieworders, CAC 190123
vieworder = obj.views.vieworder;
N_altbach = obj.views.NframesDynamic;
Nviews_header = obj.views.Nviews_header;
Nviews_header_iso = obj.views.Nviews_header_iso;
Nviews_header_seg = obj.views.Nviews_header_seg;
Nviews_header_rem = obj.views.Nviews_header_rem;
Nviews_header_nav = obj.views.Nviews_header_nav;

% djswift fmri dynamic, a hack but will work...CAC 190603 ***
recur = obj.pars.recur;
if obj.pars.vieworder_num == 96
    N_altbach = recur;
    Nviews_header_seg = Nviews_header_iso;
end

% kspaceRad sizes
size_kspaceRad = size( obj.kspaceRad);
Np = size_kspaceRad(1); Nv = size_kspaceRad(2); Ns = size_kspaceRad(3);
size_kspace_i = Np*Nv*Ns; Nviews_kspace = Nv*Ns*recur;
s_ksRad = size_kspaceRad;
kR_Npnts = s_ksRad(1)*s_ksRad(2)*s_ksRad(3);

% fid now (readout, segments, volumes, channels)
size_fid = size( obj.fid);
ss_fid = ndims( obj.fid);
% should be 4 dims if multi-channel
if ss_fid < 3    
    size_fid(3) = 1;    % single volume data
    size_fid(4) = 1;    % single channel data
elseif ss_fid < 4
    size_fid(4) = 1; 
end
Nfid = size_fid(1); Nviews = size_fid(2)*size_fid(3); Nchannels = size_fid(4);
size_fid_i = Nfid*Nviews;   % for GE data

if recur == 1 % hack, CAC 190603 ***
    if ( Nviews == Nviews_kspace )
        % everything ok
    else
        fprintf('\n'); warning( 'Data views %d not equal to expected %d views.', Nviews, Nviews_kspace);
    end
end

nrcvrs = obj.nrcvrs;
if ( nrcvrs == Nchannels )
    % everything okDscale = round( sqrt(N_altbach));
else
    fprintf('\n'); warning( 'kdata Nchannels %d not equal to nrcvrs %d ', Nchannels, nrcvrs);
end

% dynamic/Altbach information
is_altbach = obj.FLAGS.DYNAMIC;
if ~is_altbach
    N_altbach_loop = 1;
    Nviews_dynamic = Nviews_header_iso;
else
    N_altbach_loop = N_altbach;
    Nviews_dynamic = Nviews_header_seg;
end

% origin
% not yet implemented, get the orientation correct first? CAC 180927

% local param variables
red_matrix = obj.param.red_matrix;
img_matrix = obj.param.img_matrix;
nv = obj.param.nv;
nIter = obj.param.nIter;
w = obj.param.w;
nImg = obj.param.img_matrix;
OS = obj.param.oversampling;
KW = obj.param.kernelwidth;
FS = obj.param.freq_scale;
KS = obj.param.kernelsampling;
Dscale = obj.param.Dscale;

% voxel sizes for nifti, mm, sec
lro = obj.pars.lro;
vox = Dscale*lro*10/(2*Np);
total_time = obj.pars.sctime;
vox_t = total_time*60/N_altbach_loop;    % not true for view shared CAC 180927

% prepare for window
dcf = reshape( obj.dcf, Nfid, Nviews_header_iso);

%if ( DEBUG_FLAG >= 3 ); param_kdata = obj.param, end

% cGrid preparationn
% calculate kernel
kernel = Calculate_kernel( obj.param);
% details for cgrid
%nImg = obj.param.img_matrix/obj.param.freq_scale(1); %*** 121017 CAC
%nImg = obj.param.img_matrix;
cgrid_opts.kernelradius=KW;
cgrid_opts.nx=OS*nImg;
cgrid_opts.ny=OS*nImg;
cgrid_opts.nz=OS*nImg;
% details for cInvgridNew
cInvgrid_opts.kernelradius=KW;
cInvgrid_opts.nx=OS*nImg;
cInvgrid_opts.ny=OS*nImg;
cInvgrid_opts.nz=OS*nImg;

%% Dynamic Loop
obj.kspaceRad = reshape( obj.kspaceRad, Nfid, Nviews_header_iso, 3);   % setting up for dynamic
obj.fid = reshape( obj.fid, Nfid, Nviews_header_iso*recur, Nchannels); % semi flatten as needed by dynamic loop
obj.dcf = reshape( obj.dcf, Nfid, Nviews_header_iso); % semi flatten as needed by dynamic loop
obj.kspCart = complex( zeros( OS*nImg, OS*nImg, OS*nImg, N_altbach_loop, Nchannels, 'single'), zeros( OS*nImg, OS*nImg, OS*nImg, N_altbach_loop, Nchannels, 'single')); % force single, CAC 200807 ***
obj.kspAPC = complex( zeros( OS*nImg, OS*nImg, OS*nImg, N_altbach_loop, 'single'), zeros( OS*nImg, OS*nImg, OS*nImg, N_altbach_loop, 'single'));
for i_dynamic = 1:N_altbach_loop
    i_start = 1;    % defaults to 1
    i_end = Nviews_header_iso;
    if N_altbach_loop > 1
        i_start = (i_dynamic - 1)*Nviews_dynamic +1;
        i_end = i_start + Nviews_dynamic - 1;
    end
    
    % handle recur
    i_start_ks = mod( i_start, Nviews_header_iso);
    i_end_ks = mod( i_end, Nviews_header_iso);
    if i_end_ks == 0
        i_end_ks = Nviews_header_iso;
    end
    
    kspace_i = obj.kspaceRad( :, i_start_ks:i_end_ks, :);
    kspace_i = reshape( kspace_i, Nfid*Nviews_dynamic, 3);  % flatten for gridding
    
    %   check size of input_data
    if ( size( kspace_i, 1) ~=3 )
        if( size( kspace_i, 2) ~=3 )
            fprintf('\n'); error( 'In Calculate_DCF_ModPipe.m: wrong matrix size for kspace_i  ');
        end
    end
    % realign to size as required by cgrid.c
    if ( size( kspace_i, 1) == 3 )
        kspace_i = permute( kspace_i, [2 1]);
    end
    
    % gridding indices
    kspace_i = kspace_i*red_matrix*OS*FS(1);   % do by components for anisotropic? CAC 181204 ***
    sz_kspace_i = size( kspace_i, 1);
    
    %% channel loop
    for i_channel = 1:Nchannels
        fid_i = squeeze( obj.fid(:, i_start:i_end, i_channel));
        if obj.FLAGS.DMF == 1   % need to fix for %% dynamic density weight calculation navigators, 190131 CAC ***
            kspaceL_i = obj.kspaceRadLambda(:, i_start_ks:i_end_ks);    % Lambda for DMF
            fid_i =  fid_i .* kspaceL_i;
            %fid_i = angle( fid_i) .* kspaceL_i;
            %fid_i = exp( i*pi*fid_i);
        end
        
        fid_i = reshape( fid_i, Nviews_dynamic*Nfid, 1);
        
        % dynamic density calculation removed
        %   better to use windowed global density?
        
        % minimal cost to index the pre-calculated density
        % index the global density
        if ( DEBUG_FLAG >= 3 ); fprintf( '\nUsing global DCF indices %d to %d for dynamic frame %d/%d...  ', i_start_ks, i_end_ks, i_dynamic, N_altbach_loop); end
        % for normalization
        dPnorm = Nviews_header_iso/Nviews_dynamic;
        % index global iso density
        dcf_d = dPnorm * obj.dcf( :, i_start_ks:i_end_ks);
        
        if obj.switches.window_density >= 1 % window density weights
            N_hann = img_matrix - 2;
            window_for_density = hann( N_hann, 'periodic'); % Hanning window
            if obj.switches.window_density >= 2; window_for_density = sqrt( window_for_density); end % sqrt Hanning window
            if obj.switches.window_density >= 3; window_for_density = sqrt( window_for_density); end % another sqrt
            if obj.switches.window_density >= 4; window_for_density = sqrt( window_for_density); end % another sqrt
            window_for_density = window_for_density((N_hann/2+1):end);
            %window_for_density = window_for_density/window_for_density(1); % normalize?, 210501 CAC ***
            s_pad = (red_matrix - N_hann/2);
            window_for_density = padarray(window_for_density, s_pad, 0, 'post');
            if ( DEBUG_FLAG >= 3 ); fprintf( '\nwindowing density weights with window #%d\n', obj.switches.window_density); end
            if ( DEBUG_FLAG >= 4 ); figure( 'Name', 'window_for_density'); plot( window_for_density); end
            dcf_d = repmat( window_for_density, 1, Nviews_dynamic) .* dcf_d;
        end
        
        % flatten dP
        dcf_d = reshape( dcf_d, Nfid*Nviews_dynamic, 1);
        
        %% Grid
        cur_size_fid_i = size( fid_i, 1);
        if ( DEBUG_FLAG >= 2 ); fprintf( 'Starting forward grid of %d k-space points for dynamic frame %d/%d channel %d/%d, ', cur_size_fid_i, i_dynamic, N_altbach_loop, i_channel, Nchannels); end
        if ( DEBUG_FLAG >= 2  & ( i_channel < Nchannels | i_dynamic < N_altbach_loop) );  fprintf( '\n'); end
        %kR_norm = (1.0/1.14)/OS^3;  % 1.14 is amplitude with OS = 2.0, KW = 4.0
        %kR_norm = (1.00/2.384896)/OS^3;    % 2.41 is amplitude with OS = 1.25, KW = 2.5
        kR_norm = 1.00/(2.391*OS^3); % 2.391 is amplitude with OS = 1.25, KW = 2.5, average at slice 65 vienum 128, CAC 210501 ***
        %kR_norm = 1.00;    % for kR_norm meausrement
        fid_i = kR_norm * fid_i .* dcf_d;
        
        obj.kspCart(:, :, :, i_dynamic, i_channel) = gridding3D_mex( fid_i, kspace_i, kernel, cgrid_opts);
        
        %  * grid_volume = grid3_MAT(0:data, 1:crds, 2:weights, 3:effMkern,kbu] = calckbkernel(kwidth,overgridfactor,klength)tx, 4:numThreads)
        %  *
        %  * REQUIRED:
        %  *  data: N-D double-precision matrix >=2D, fastest varying dimension is length 2
        %  *
        %  *  coords: N-D double-precision matrix >=2D, fastest varying dimension is length 3
        %  *              trajectory coordinate points scaled between -0.5 to 0.5
        %  *
        %  *  weights: N-D double-precision matrix >=1D, size(coords)/3 == size(weights) == size(data)/2
        %  *          if no weighting is desired then create a matrix of the appropriate size and set values
        %  *          to unity.
        %  *
        %  *  effMtx: (integer) the length of one side of the grid matrix, range >=1
        %  *
        %  *  numThreads: (integer) 1 or 8, the code will run either non-threaded or threading on grid octants
        %  */
        
        
    end
    
    % calculate apodization correction
    % move to doDCF and cache with dcf!!! *** CAC 190302
    if obj.switches.calc_apc
        fid_apc = complex( zeros( size( fid_i), 'single')); fid_apc(1) = complex(1.0);
        obj.kspAPC(:, :, :, i_dynamic) = gridding3D_mex( fid_apc, kspace_i, kernel, cgrid_opts);
    end
end

size_ksc = size( obj.kspCart);
if ( DEBUG_FLAG >= 3 )
    fprintf( '\n%s DEBUG kspCart array size: %d %d %d %d %d\n', mfilename, size_ksc')
end

size_ksa = size( obj.kspAPC);
if ( DEBUG_FLAG >= 3 )
    fprintf( '\n%s DEBUG kspAPC array size: %d %d %d %d\n', mfilename, size_ksa')
end

%plot( abs( squeeze( obj.kspCart(end/2, end/2, :))))

obj.kspaceRad = reshape( obj.kspaceRad, Nfid, Nv, [], 3); % return to original shape
obj.dcf = reshape( obj.dcf, Nfid, Nv, []); % return to original shape
obj.fid = reshape( obj.fid, Nfid, Nv, [], Nchannels); % return to original shape

% clear out memory
%what can be cleared? (won't help anyways? *** CAC 210405)

% elapsed time
if DEBUG_FLAG >= 2
    %fprintf('kdata array size: %d %d %d %d\n', size(obj.kdata));
    toc( kdata_timing);
end

return;
