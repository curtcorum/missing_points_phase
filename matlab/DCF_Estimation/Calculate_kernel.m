function [kernel] = Calculate_kernel(PARAM)

if(~isfield(PARAM,'oversampling'))
    error('In gridding3D_mex: PARAM missing field - oversampling');
end
if(~isfield(PARAM,'red_matrix'))
    error('In gridding3D_mex: PARAM missing field - red_matrix');
end
if(~isfield(PARAM,'kernelwidth'))
    error('In gridding3D_mex: PARAM missing field - kernelwidth');
end
if(~isfield(PARAM,'kernelsampling'))
    error('In gridding3D_mex: PARAM missing field - kernelsampling');
end

% define kernel 
nsamp_pos=ceil(PARAM.red_matrix*PARAM.oversampling);  % number of samples from the center out;
readout_spacing=.5/(PARAM.oversampling*PARAM.red_matrix);
nsmap_pos2=nsamp_pos+ceil(PARAM.oversampling*PARAM.kernelwidth);
finalgrid=(-nsmap_pos2:1:nsmap_pos2)*readout_spacing;

nsamp_pos=PARAM.red_matrix; %floor(PARAM.matrix/2);
nsamp_neg=PARAM.red_matrix*2-nsamp_pos; % number of samples over the center ;
nsamp=max([nsamp_pos nsamp_neg]);
[kernel,kernelgrid,final_short]=kernel_setup_1D_x2(finalgrid,PARAM.kernelwidth,nsamp,...
    PARAM.oversampling ,PARAM.kernelsampling);
kernel=real(kernel);
end