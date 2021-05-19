function [dP, d0] = Calculate_DCF_ModPipe_gpuNUFFT( kspaceRadt, param, nIter)
% function [dP, d0] = Calculate_DCF_ModPipe_gpuNUFFT( kspaceRadt, param, nIter)
%   kspaceRadt  the k-space co-ordeinates for the sampling scheme (transposed so gradients first)
%   param       structure containing the gridding parameters
%   nIter       number of iterations, will use value in param if nIter = -1
%
%   dP          density correction
%   d0          initial density correction estimate (one iteration)
%
% Curt Corum
% Champaign Imaging LLC 3/1/2019
%
% converted to use gpuNUFFT

% use cGrid ModPipe3D if gpuNUFFT not available? *** CAC 190301

% check size of input_data
if( size( kspaceRadt, 1) ~= 3 )
    if( size( input_data, 2) ~=3 )
        error( 'wrong matrix size for input_data, not 3xN')
    end
end

% calculate kernel 
%kernel = Calculate_kernel(PARAM);
% scaling to make area under the curve 1
%scale = 1/sum(kron(kernel,kron(kernel,kernel)));
%kernel = kernel*scale;

%gpuNUFFT wants opposite of cGrid, checks inside gpuNUFFT CAC 190301
% realign to size as required by cgrid.c
% if(size(input_data,1)==3)
%     input_data=permute(input_data,[2 1]);
% end

% local param variables
red_matrix = param.red_matrix;
img_matrix = param.img_matrix;
nv = param.nv;
nIter = param.nIter;
w = param.w;
nImg = param.img_matrix;
OS = param.oversampling;
KW = param.kernelwidth;
FS = param.freq_scale;
KS = param.kernelsampling;

% local variables
SW = 8;
imageDim = OS*nImg*[1 1 1];


% details for cgrid
%nImg = PARAM.img_matrix*PARAM.freq_scale(1);
% nImg = PARAM.img_matrix;
% cgrid_opts.kernelradius=PARAM.kernelwidth;
% cgrid_opts.nx=PARAM.oversampling*nImg;
% cgrid_opts.ny=PARAM.oversampling*nImg;
% cgrid_opts.nz=PARAM.oversampling*nImg;

% details for cInvgridNew
% cInvgrid_opts.kernelradius=PARAM.kernelwidth;
% cInvgrid_opts.nx=PARAM.oversampling*nImg;
% cInvgrid_opts.ny=PARAM.oversampling*nImg;
% cInvgrid_opts.nz=PARAM.oversampling*nImg;

% gridding indices
%kmap = input_data*PARAM.red_matrix*PARAM.oversampling*PARAM.freq_scale(1);
% gpuNUFFT wants k-space coordinates in +-0.5
if ( max( abs( kspaceRadt), [], 'all') > 0.5 )
    warning( 'kspaceRadt scaling not +-0.5');
end

% function calls for H(regrid: Cartesian - Radial) and Ht(grid: Radial - Cartesian)
% fnHt = @(ksp) gridding3D_mex(ksp,kmap,kernel,cgrid_opts);
% fnH = @(ksp) regrid3D_mex(ksp,kmap,kernel,cInvgrid_opts);

        % function FT = gpuNUFFT(k,w,osf,wg,sw,imageDim,sens,varargin)
        %
        %     k -- k-trajectory, scaled -0.5 to 0.5
        %          dims: 2 (3) ... x y (z)
        %                N ... # sample points
        %     w -- k-space weighting, density compensation
        %     osf -- oversampling factor (usually between 1 and 2)
        %     wg -- interpolation kernel width (usually 3 to 7)
        %     sw -- sector width to use
        %     imageDim -- image dimensions [n n n]
        %     sens -- coil sensitivity data
        %     varargin
        %        opt  -- true/false for atomic operation (default true)
        %             -- true/false for using textures on gpu (default true)
        %             -- true/false for balanced operation (default true)
        %
        %  FT -- gpuNUFFT operator
        %
        %  A. Schwarzl, Graz University of Technology
        %  F. Knoll, NYU School of Medicine
        % 
        % Now need to pass w.^2 for density if doing grid once!
        % see:  https://github.com/andyschwarzl/gpuNUFFT/issues/59#issuecomment-320018001

if size( kspaceRadt, 1) ~= 3; error( ' kspaceRadt size is incorrect'); end
m = size( kspaceRadt, 2);
w = ones( m, 1);
FT = gpuNUFFT( kspaceRadt, w, OS, KW, SW, imageDim, []);

% progress
%fprintf( '\nStarting Pipe/Menon DCF with %d iterations for %d k-space points...  ', nIter, size( input_data, 1));

% call ModPipe3D
[dP, d0] = ModPipe3D_gpuNUFFT( FT, FT', m, nIter, []);

end



