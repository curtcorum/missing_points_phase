function kdata = cubeObject( kspace, varargin)
%function kdata = cubeObject( kspace, gain, kmax, [Lx Ly Lz], [Dx Dy Dz], R, Df, TA, T2)
% generate k-space data point(s) for a cube object
%
%   kdata       output array (N) of data values for each k-space coordinate 
%
%   kspace      array (N*3) of k-space points
%   gain        amplitude of object
%   kmax        max kspace radius (abs)
%   [Lx Ly Lz]  length of object (FOV units)
%   [Dx Dy Dz]  displacement of object (FOV units)
%   R           rotation matrix of object
%   Df          off resonance parameter, Hz
%   TA          aquisition time (center out in s)
%   T2          T2 in s
%
%   defaults: gain = 1; L = [.5 .5 .5]; D = [0 0 0]; R = eye(3);
%
% Curt Corum, Champaign Imaging LLC, 5/15/2019
%   200809 CAC - off resonance and T2 decay added
%   210224 CAC - displacement in unrotated coordinates

% defaults:
gain = 1; kmax = 1; L = [.5 .5 .5]; D = [0 0 0]; R = eye(3);
Df = 100; TA = 4.096e-3; T2 = 1.0e-3;

% varargin handling
if nargin > 1; gain = varargin{1}; end
if nargin > 2; kmax = varargin{2}; end
if nargin > 3; L = varargin{3}; end
if nargin > 4; D = varargin{4}; end
if nargin > 5; R = varargin{5}; end
if nargin > 6; Df = varargin{6}; end
if nargin > 7; TA = varargin{7}; end
if nargin > 8; T2 = varargin{8}; end
if nargin > 9; error('too many arguments'); end

% input checking
sz_kspace = size( kspace);
n_points = prod( sz_kspace(1:(end-1)));
if sz_kspace( end) ~= 3
    error( 'kspace must be an array of 3-vectors')
end

% flattening and preparing kspace
kspace = reshape( kspace, n_points, 3);
kspace = kspace * R;

% displacements in unrotated coordinates
D = D * R;

%argt = 1/(norm( [ 1 1 1])) * kspace * [1 1 1]'; argt = abs( argt);
argt = sqrt( kspace(:, 1) .* kspace(:, 1) + kspace(:, 2) .* kspace(:, 2) + kspace(:, 3) .* kspace(:, 3));
argx = .5 * kmax * kspace(:, 1) * L(1);
argy = .5 * kmax * kspace(:, 2) * L(2);
argz = .5 * kmax * kspace(:, 3) * L(3);


%size(D')
%size(kspace)
%size( kspace * D')

% calculation
% if kspace is single then kdata will be single
kdata = (gain * L(1)*sinc( argx)) .* (L(2)*sinc( argy)) .* (L(3)*sinc( argz)) .* exp( i*2*pi*(.5*kmax*kspace*D') - argt*TA*(1/T2-i*2*pi*Df));
% return shaped kdata
kdata = reshape( kdata, [sz_kspace(1:(end-1)), 1]);

return
