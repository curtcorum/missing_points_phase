function [y,kernelgrid,final_short]=kernel_setup_1D(grid_spacing,kernelwidth,nsamp,oversampling,KB_kernel_res)
%grid_spacing=finalgrid;kernelwidth=PARAM.kernelwidth;nsamp=.5./PARAM.red_m
%atrix;oversampling= PARAM.oversampling;KB_kernel_res=PARAM.kernelsampling

     KB_sampling=round(KB_kernel_res*kernelwidth*oversampling)+1;

LIST=grid_spacing; 
final_data=zeros(size(LIST,2),1);

datum=[0];%sampling(10);

[A]=find(abs((LIST)-(datum(1)))==min(abs((LIST)-(datum(1)))));A=A(1);
deltax=-abs((LIST(A))-(datum(1)))*nsamp*oversampling;
A=ceil((size(LIST,2)+1)/2);   %% NOTE the ceil
ROI_lb_x=A-ceil(oversampling*kernelwidth/2);
     ROI_ub_x=A+ceil(oversampling*kernelwidth/2);

readout_int=linspace(LIST(ROI_lb_x),LIST(ROI_ub_x),KB_sampling);

Xshort=readout_int;

X=LIST(ROI_lb_x:ROI_ub_x);

final_short.X=X;

u=Xshort*1;

%%%%%%%%%%%%
a = oversampling;
w = kernelwidth;	
beta = pi*sqrt( w^2/a^2*(a-0.5)^2-0.8 );	% From Beatty et al.
beta=beta;
x= beta*sqrt(1-(2*u*nsamp/w).^2);	% Argument - see Jackson '91.

y= (besseli(0,x.^2)./w);% the saptial variable is x^2, instead of linear in x
%size(y)
y=y./y( (size(y,1)+1)/2,  (size(y,2)+1)/2,  (size(y,3)+1)/2);  % normalization

kernelgrid.X=linspace(1,size(final_short.X,2),KB_sampling);



y=real(y(ceil(end/2):end));
kernelgrid.X=linspace(1,(size(final_short.X,2)-1)/2+1,size(y,2));
%kernelgrid.X=linspace(1,a+1,size(y,2));

