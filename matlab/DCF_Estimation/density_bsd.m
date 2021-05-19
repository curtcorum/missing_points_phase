function [d, d0, w] = density_bsd( fnH, fnHt, n, m, w, iter)
%function [d, d0] = density_bsd( fnH, fnHt, n, m, w, iter)
%   fnH is the interpolation function. fnHt is the adjunct funtion.
%   n is size of on dimension of oversampled cartesian matrix (os*Nx)
%   m is number (total points) of non-cartesian samples
%   w is a scalar controlling the constraint. -1 calculates automatic estimate.
%   iter is the number of iterations.
%   Suggested values are w = 2*max(max(H)) and iter = 10. (for 2d)
%   
% originally from: Bydder M, Samsonov AA, Du J. Evaluation of optimal density weighting for regridding. Magnetic resonance imaging. 2007 Jun 1;25(5):695-702, Appendix A.
% http://dx.doi.org/10.1016/j.mri.2006.09.021
%
% Curt Corum, Champaign Imaging LLC, 1/13/2019
%
% CAC 190113    Modified to use RV DCF_Estination function calls for H and Ht

if w == -1
    w = 2*max( max( fnH( ones( n, n, n))));
end

d = 1./(fnH( fnHt( ones( m, 1))));
d0 = d; % save grid once density

g = 1./d+w^2;
r = fnH( ones( n, n, n) - fnHt( d));
s = r./g;
new = r'*s;
for i= 1:iter
    q = fnH( fnHt( s))+w^2*s;
    alpha = new/(s'*q);
    d = d+alpha*s;
    r = r-alpha*q;
    q = r./g;
    old = new;
    new = r'*q;
    s = q+(new/old)*s;
end

end