function obj = doImages( obj)
% ngfnRecon.doImages
% takes kspCart and produces images
%   C. Corum 2/26/2019
%
%   CAC 190226 edited out of doKdata
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

if ( DEBUG_FLAG >= 2 ); img_timing = tic; end     % start timer

%% Switches

%debug
%switches = obj.switches

%% Set up variables
% needs a lot of cleaning up in here, doKdata and doDCF. *** CAC 190226

% paths
startPath = obj.startPath;
savePath = obj.savePath;

% sizes
size_ksc = size( obj.kspCart);
if ( ndims( obj.kspCart) < 4 ); size_ksc(4) = 1; end
if ( ndims( obj.kspCart) < 5 ); size_ksc(5) = 1; end
size_ksa = size( obj.kspAPC);
if ( ndims( obj.kspAPC) < 4 ); size_ksa(4) = 1; end
N_altbach = size_ksc(4);
Nchannels = size_ksc(5);

% Don't forget about Navgators coordinate transforms..., *** CAC 190206

% dynamic/Altbach information
is_altbach = obj.FLAGS.DYNAMIC;
if ~is_altbach
    N_altbach_loop = 1;
else
    N_altbach_loop = N_altbach;
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
Np = obj.param.red_matrix;
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

% crop
if obj.switches.crop == 1
    cropStart = floor( (nImg * obj.param.oversampling - nImg)/2) + 1;
    cropEnd   = floor( (nImg * obj.param.oversampling + nImg)/2);
    nImg_final = nImg;
elseif obj.switches.crop == 2
    cropStart = floor( (nImg * obj.param.oversampling - nImg/2)/2) + 1;
    cropEnd   = floor( (nImg * obj.param.oversampling + nImg/2)/2);
    nImg_final = nImg/2;
else
    %imageOutShiftCrop = imageOutShift;
    cropStart = 1;
    cropEnd = nImg * obj.param.oversampling;
    nImg_final = cropEnd;
end

%if ( DEBUG_FLAG >= 3 ); param_img = obj.param, end

%% Dynamic Loop
% initialize img arrays
% do fft etc. in place with kspCart to save memory? *** CAC 210405
obj.img = complex( zeros( nImg_final, nImg_final, nImg_final, N_altbach_loop, Nchannels, 'single'), zeros( nImg_final, nImg_final, nImg_final, N_altbach_loop, Nchannels, 'single'));
obj.img_sos = complex( zeros( nImg_final, nImg_final, nImg_final, N_altbach_loop, 'single'), zeros( nImg_final, nImg_final, nImg_final, N_altbach_loop, 'single'));
for i_dynamic = 1:N_altbach_loop
    
    % channel loop
    for i_channel = 1:Nchannels
        % fft
        if ( DEBUG_FLAG >= 2  & ( i_channel <= Nchannels | i_dynamic < N_altbach_loop ) & obj.switches.do_fft & ( i_channel ~= 1 ));  fprintf( '\n'); end
        if ( DEBUG_FLAG >= 2 ) & obj.switches.do_fft; fprintf( '3dFFT of size %d^3 for dynamic frame %d/%d channel %d/%d, ', nImg * obj.param.oversampling, i_dynamic, N_altbach_loop, i_channel, Nchannels); end
        if obj.switches.fft_shift; kspCartShift = fftshift( obj.kspCart(:, :, :, i_dynamic, i_channel)); else kspCartShift = obj.kspCart(:, :, :, i_dynamic, i_channel); end
        if obj.switches.do_fft; img = fftn( kspCartShift); else img = kspCartShift; end
        clear kspCartShift;
        if obj.switches.fft_shift_2; img = fftshift( img); end        

        % apodization correction            *** Calculation can be done once??? *** CAC 210326
        if obj.switches.apc
            %fprintf( '...FFT and Shifting...\n');
            kspAPC = fftshift( obj.kspAPC(:, :, :, i_dynamic));
            apcOut = fftn( kspAPC);
            % normalize to center of fov so correponds to uncorrected, *** CAC 170309 ***
            apcOutNorm = apcOut(1, 1, 1);
            if ( DEBUG_FLAG >= 3 ); fprintf( '\nAPC normalization: %g...', apcOutNorm); end
            apcOut = apcOut/apcOutNorm;
            apcOut = fftshift( apcOut);
            apcOut(~(abs( apcOut))) = complex( 1.0);
            img = img ./ apcOut;
        end        
        
        img = img(cropStart:cropEnd, cropStart:cropEnd, cropStart:cropEnd);
        
        % if switches.save_to_obj % save memory
        obj.img(:, :, :, i_dynamic, i_channel) = img;
        % end
        
        % debug display of AvgAbs
        %imageOutAvgAbs = mean( mean( mean( abs(imageOut))));
        % print message and add backspaces to Nbsp
        %fprintf(repmat('\b',1,Nbsp));
        %msg = sprintf( 'Average value of abs pixels %g', imageOutAvgAbs);
        %Nbsp = length(msg);
        %fprintf('%s', msg);
        
        % save cropped image as img files in nifti.d
        if obj.switches.nifti
            magFile = 'mag.img';
            phaseFile = 'phase.img';
            dirName = 'nifti.d'; cd( savePath);
            if ( DEBUG_FLAG >= 2 ); fprintf( 'Saving %s and %s files in %s, ', magFile, phaseFile, dirName); end
            Write_Raw_Sgl( magFile, phaseFile, dirName, img);  % Write_Raw_Sgl creates dirName
        end
        
        % save coil sensitivities in nifti_csm.d
        if obj.switches.nifti_csm & (i_dynamic == 1)
            magFile = 'mag.img';
            phaseFile = 'phase.img';
            dirName = 'nifti_csm.d'; cd( savePath);
            if ( DEBUG_FLAG >= 2 ); fprintf( 'Saving %s and %s files in %s, ', magFile, phaseFile, dirName); end
            Write_Raw_Sgl( magFile, phaseFile, dirName, obj.img_csm(:, :, :, i_channel));  % Write_Raw_Sgl creates dirName
        end     
        
        % accumulate sos
        %size_i = size( img)
        %size_is = size( obj.img_sos)
        obj.img_sos(:, :, :, i_dynamic) = obj.img_sos(:, :, :, i_dynamic) + abs( img .* conj( img)); % obj.img_sos = img_sos, *** CAC 190919 do this at some point?
        
    end     % end channel loop
    
    %save img_rsos in nifti_sos.d
    if obj.switches.nifti_sos
        magFile = 'mag.img';
        phaseFile = 'NONE';
        dirName = 'nifti_sos.d'; cd( savePath);
        if ( DEBUG_FLAG >= 2 ); fprintf( '\nSaving %s file in %s, ', magFile, dirName); end
        obj.img_sos(:, :, :, i_dynamic) = sqrt( obj.img_sos(:, :, :, i_dynamic));
        Write_Raw_Sgl( magFile, phaseFile, dirName, obj.img_sos(:, :, :, i_dynamic));  % Write_Raw_Sgl creates dirName
    end
    
end     % end dynamic (Altbach) loop


%% Final save and exit
% Need to autosave relavant nifti header information for fMRI and reference
% images CAC 180927

% save nifti hdr files in nifti.d
if obj.switches.nifti
    if ( DEBUG_FLAG >= 2 ); fprintf( 'nifti header file, '); end
    
    if is_altbach
        imgsize = [nImg_final nImg_final nImg_final N_altbach_loop Nchannels];
        voxel_size = [vox vox vox vox_t 1];
    else
        imgsize = [nImg_final nImg_final nImg_final Nchannels];
        voxel_size = [vox vox vox 1];
    end
    origin = [0 0 0];       %read from Pfile?
    datatype = 16;	% float32
    description = strcat('ngfn_', version);  % header version
    maxval_mag = max( abs( obj.img), [], 'all');	minval_mag = 0.00;
    hdr_mag = make_nii_hdr(imgsize, voxel_size, origin, datatype, description, maxval_mag, minval_mag);
    hdr_mag.dime.xyzt_units = 2 + 8;    % mm and sec
    maxval_phase = max( angle( obj.img), [], 'all');	minval_phase = min( angle( obj.img), [], 'all');
    hdr_phase = make_nii_hdr( imgsize, voxel_size, origin, datatype, description, maxval_phase, minval_phase);
    hdr_phase.dime.xyzt_units = 2 + 8;    % mm and sec
    dirName = 'nifti.d'; cd( savePath); cd( dirName); % already exists from img save
    magHdr = 'mag.hdr';	phaseHdr = 'phase.hdr';
    % open the files and write nii headers
    magHdrFid = fopen( magHdr, 'a+', 'ieee-be'); phaseHdrFid = fopen( phaseHdr, 'a+', 'ieee-be');
    save_nii_hdr( hdr_mag, magHdrFid); save_nii_hdr( hdr_phase, phaseHdrFid);
    %close the files
    stDebugr = fclose( magHdrFid); stDebugi = fclose( phaseHdrFid);
    obj.hdr = hdr_mag;
    cd( startPath);
end

% Save nii and nii.gz testing, *** CAC 190708
%dirName = 'nifti.d'; cd( savePath); cd( dirName);
%save_compressed_nii(nii_struct, output_path, untouch)

if obj.switches.nifti_sos
    if ( DEBUG_FLAG >= 2 ); fprintf( 'nifti_sos header file, '); end
    
    if is_altbach
        imgsize = [nImg_final nImg_final nImg_final N_altbach_loop];
        voxel_size = [vox vox vox vox_t];
    else
        imgsize = [nImg_final nImg_final nImg_final];
        voxel_size = [vox vox vox];
    end
    origin = [0 0 0];       %read from Pfile?
    datatype = 16;	% float32
    description = strcat('ngfn_', version);  % header version
    maxval_mag = max( abs( obj.img_sos), [], 'all');	minval_mag = 0.00;
    hdr_mag = make_nii_hdr(imgsize, voxel_size, origin, datatype, description, maxval_mag, minval_mag);
    hdr_mag.dime.xyzt_units = 2 + 8;    % mm and sec
    %maxval_phase = max( angle( obj.img), [], 'all');	minval_phase = min( angle( obj.img), [], 'all');
    %hdr_phase = make_nii_hdr( imgsize, voxel_size, origin, datatype, description, maxval_phase, minval_phase);
    dirName = 'nifti_sos.d'; cd( savePath); cd( dirName); % already exists from img save
    magHdr = 'mag.hdr';	%phaseHdr = 'phase.hdr';
    % open the files and write nii headers
    magHdrFid = fopen( magHdr, 'a+', 'ieee-be'); %phaseHdrFid = fopen( phaseHdr, 'a+', 'ieee-be');
    save_nii_hdr( hdr_mag, magHdrFid); %save_nii_hdr( hdr_phase, phaseHdrFid);
    %close the files
    stDebugr = fclose( magHdrFid); %stDebugi = fclose( phaseHdrFid);
    obj.hdr_sos = hdr_mag; % need to save the dynamic series in memory, make movie, *** CAC 190226
    cd( startPath);
end

% save nifti hdr files in nifti_csm.d
if obj.switches.nifti_csm
    if ( DEBUG_FLAG >= 2 ); fprintf( 'nifti_csm header file, '); end
    imgsize = [nImg nImg nImg Nchannels];
    voxel_size = [vox vox vox 1];
    origin = [0 0 0];       %read from Pfile?
    datatype = 16;	% float32
    description = strcat('ngfn_', version);  % header version
    maxval_mag = max( abs( obj.img_csm), [], 'all');	minval_mag = 0.00;
    hdr_mag = make_nii_hdr(imgsize, voxel_size, origin, datatype, description, maxval_mag, minval_mag);
    hdr_mag.dime.xyzt_units = 2 + 8;    % mm and sec
    maxval_phase = max( angle( obj.img_csm), [], 'all');	minval_phase = min( angle( obj.img_csm), [], 'all');
    hdr_phase = make_nii_hdr( imgsize, voxel_size, origin, datatype, description, maxval_phase, minval_phase);
    hdr_phase.dime.xyzt_units = 2 + 8;    % mm and sec
    dirName = 'nifti_csm.d'; cd( savePath); cd( dirName); % already exists from img save
    magHdr = 'mag.hdr';	phaseHdr = 'phase.hdr';
    % open the files and write nii headers
    magHdrFid = fopen( magHdr, 'a+', 'ieee-be'); phaseHdrFid = fopen( phaseHdr, 'a+', 'ieee-be');
    save_nii_hdr( hdr_mag, magHdrFid); save_nii_hdr( hdr_phase, phaseHdrFid);
    %close the files
    stDebugr = fclose( magHdrFid); stDebugi = fclose( phaseHdrFid);
    obj.hdr_csm = hdr_mag;
    cd( startPath);
end

cd( startPath);
if DEBUG_FLAG >= 2 || obj.switches.png == 1
    figure( 'name', obj.fidName);
    hold on; title( 'doKdata DEBUG center z slice magnitude'); colormap( 'Gray'); colorbar;
    Image_preview = squeeze( obj.img_sos(:, :, 1+(end/2), 1));
    imagesc( Image_preview);
    hold off;
    maxval_mag = max( max( Image_preview));
    Image_preview = 1.05* Image_preview/maxval_mag; %scale to 95% of max
    Image_preview = flip( Image_preview, 1); %correct the (y) orientation relative to screen and nifti
    Software_string = sprintf( 'ngfnRecon version = %s', version);
    cd( savePath);
    %fprintf( '\npreview.png...');
    imwrite( Image_preview, 'preview.png', 'PNG', 'BitDepth', 16, 'InterlaceType', 'adam7', ...
        'ResolutionUnit', 'meter', 'XResolution', 1e3/vox, 'YResolution', 1e3/vox, 'SignificantBits', 12, ...
        'Source',  savePath, 'Copyright', 'Champaign Imaging LLC, all rights reserved', 'Software', Software_string);
    cd( startPath);
end

% elapsed time
if DEBUG_FLAG >= 2
    %fprintf('kdata array size: %d %d %d %d\n', size(obj.kdata));
    toc( img_timing);
end

return;

% if ( DEBUG_FLAG >= 2  & ( i_channel < Nchannels | i_dynamic < N_altbach_loop )) & obj.switches.do_fft;  fprintf( '\n'); end
 
