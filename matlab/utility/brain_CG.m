%function obj = doKdata_cg( obj)
% ngfnRecon.doKdata_cg
% Grids radial k-space data to cartesian k-space
%   C. Corum 5/3/2015
%
%   CAC 190226  separate DCF calculation, gridding, image manipulation and saving.

%   korig - the kspace coordinates
%
% gpuNUFFT
%   mtx_hr - kspace array size (readout)
%   ones(size(dcf(:)) - the "dummy" all ones density correction function
%   osf -- oversampling factor (usually between 1 and 2)
%   wg -- interpolation kernel width (usually 3 to 7)
%   imageDim -- image dimensions [n n n]    [1 1 1]*mtx_hr
%   sens -- coil sensitivity data   []
%   varargin 
%       opt  -- true/false for atomic operation (default true)
%            -- true/false for using textures on gpu (default true)
%            -- true/false for balanced operation (default true)


korig = 0.5*korig*mtx_hr;

FT = gpuNUFFT( transpose(korig)/mtx_hr, ones( size( dcf(:))), osf, wg, sw, [1 1 1]*mtx_hr, [], true);

mo_crpt = FT'*(kdata(:)); % data agreement in gridded k-space

ATA = @(x) (reshape( FT'*(FT*reshape( x, [mtx_hr, mtx_hr, mtx_hr])), [mtx_hr^3, 1])); % function handle

tic();
% X = pcg(A,B,TOL,MAXIT,M1,M2,X0) specifies the initial guess. If X0 is [] then pcg uses the default, an all zero vector.
[x1,flag1,rr1,iter1,rv1] = pcg( ATA, mo_crpt(:), 1e-5, 300);
toc();