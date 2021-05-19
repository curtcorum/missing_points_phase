function [dP, d0, dmax, dmin, kernel] = Calculate_DCF_ModPipe( input_data, PARAM, nIter)
% function [dP, d0, dmax, dmin] = Calculate_DCF_ModPipe( input_data, PARAM, nIter)
% dP calculated DCF
% dP0 initial DCF estimate
% dmax, dmin convergence info (largest and smallest non zero values vs iteration number)

% check size of input_data
if ( size( input_data, 1) ~= 3 )
    if ( size( input_data, 2) ~= 3 )
        error( 'wrong matrix size for input_data')
    end
end

% calculate kernel 
kernel = Calculate_kernel( PARAM);
% scaling to make area under the curve 1
% scale = 1 / sum( kron( kernel,  kron( kernel, kernel)));
% kernel = kernel * scale;

% realign to size as required by cgrid.c
if ( size( input_data, 1) == 3 )
    input_data = permute( input_data, [2 1]);
end

% details for cgrid
%nImg = PARAM.img_matrix*PARAM.freq_scale(1);
nImg = PARAM.img_matrix;
cgrid_opts.kernelradius=PARAM.kernelwidth;
cgrid_opts.nx=PARAM.oversampling*nImg;
cgrid_opts.ny=PARAM.oversampling*nImg;
cgrid_opts.nz=PARAM.oversampling*nImg;

% details for cInvgridNew
cInvgrid_opts.kernelradius=PARAM.kernelwidth;
cInvgrid_opts.nx=PARAM.oversampling*nImg;
cInvgrid_opts.ny=PARAM.oversampling*nImg;
cInvgrid_opts.nz=PARAM.oversampling*nImg;

% gridding indices
kmap = input_data*PARAM.red_matrix*PARAM.oversampling*PARAM.freq_scale(1);
% function calls for H(regrid: Cartesian - Radial) and Ht(grid: Radial - Cartesian)
fnHt = @(ksp) gridding3D_mex(ksp,kmap,kernel,cgrid_opts);
fnH = @(ksp) regrid3D_mex(ksp,kmap,kernel,cInvgrid_opts);
% progress
%fprintf( '\nStarting Pipe/Menon DCF with %d iterations for %d k-space points...  ', nIter, size( input_data, 1));
% call ModPipe3D
[dP, d0, dmax, dmin] = ModPipe3D( fnH, fnHt, size( input_data, 1), nIter, []);
end



