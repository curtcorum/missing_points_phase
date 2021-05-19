function [d, d0, dmax, dmin] = ModPipe3D( fnH, fnHt, m, iter, d0)
%function [d, d0] = ModPipe3D( fnH, fnHt, m, iter, d0)
% Iteratively calculates the sample density 'd' using the Pipe-Menon method
%  and optionally outputs 'd0' for future use.
%   fnH     function handle for regrid: Cartesian - Radial
%   fnHt    function handle for grid: Radial - Cartesian
%   m       number of k-space points
%   iter    number of iterations
%   d0      inital estimate (if available)
% Original code Rajesh V.  ~2010


% make this function a method of djswift? *** CAC 150503
%djswift.ModPipe3D

%DEBUG_FLAG = obj.FLAGS.DEBUG;
%DEBUG_FLAG = 0;
DEBUG_FLAG = 3;
%DEBUG_FLAG = 4;

% track convergence
dmax = zeros( iter+1, 1); dmin = dmax;


% Density weighting for regridding using Pipe's method
%
% ****************** Must change funtion handle call to function call for Octave ?????????? **************************
if( isempty( d0))
    d0       = fnH( fnHt( ones( m, 1)));
    d0(d0>0) = 1./d0(d0>0); % initial d
end
d = d0;
dmax(1) = max( d0(d0~=0), [], 'all'); dmin(1) = min( d0(d0~=0), [], 'all'); % debugging, CAC *** 190302
if ( DEBUG_FLAG >= 3 )
    fprintf( '\nd0max = %d d0min = %d', dmax(1), dmin(1));
end
% d0_32k = d0;
% save dcfdata.mat d0_32k -append;
% clear d0_32k;

% ****************** Must change funtion handle call to function call for Octave ?????????? **************************
for it = 1:iter
    dtmp = fnH( fnHt( d));
    condz = dtmp ~= 0;
    d(condz) = d(condz)./dtmp(condz);
    dmax(it+1) = max( d(condz), [], 'all'); dmin(it+1) = min( d(condz), [], 'all'); % debugging, CAC *** 190302
    if ( DEBUG_FLAG >= 3 )
        fprintf( '\ndmax = %d  dmin = %d iteration: %d/%d', dmax(it+1), dmin(it+1), it, iter);
    end
    
    % save intermediate dcf's
    %     if(sum(it==[1 5 10 20 50 100])~=0)
    %         eval(sprintf('dcf_%d_32k = d;',it))
    %         save dcfdata.mat -append dcf_*;
    %         clear dcf_*;
    %     end
end

if ( DEBUG_FLAG >= 3 )
    fprintf( '\n');
end

if any( d<0)
    warning( '*** WARNING - negative values in density ***')
end

end
