function kdata = gaussObject( kspace, varargin)
%function kdata = gaussObject( kspace, gain, kmax, [Lx Ly Lz], [Dx Dy Dz], R)
% generate k-space data point(s) for a cube object
%
%   kdata       output array (N) of data values for each k-space coordinate 
%
%   kspace      array (N*3) of k-space points
%   gain        amplitude of object
%   kmax        max kspace radius (abs)
%   [Lx Ly Lz]  length of object (FWHM)
%   [Dx Dy Dz]  displacement of object
%   R           rotation matrix of object (before displacement)
%
%   defaults: gain = 1; L = [.5 .5 .5]; D = [0 0 0]; R = eye(3);
%
% Curt Corum, Champaign Imaging LLC, 5/16/2019

% defaults:
gain = 1; kmax = 1; L = [.5 .5 .5]; D = [0 0 0]; R = eye(3);

% varargin handling
if nargin > 1; gain = varargin{1}; end
if nargin > 2; kmax = varargin{2}; end
if nargin > 3; L = varargin{3}; end
if nargin > 4; D = varargin{4}; end
if nargin > 5; R = varargin{5}; end
if nargin > 6; error('too many arguments'); end

% input checking
sz_kspace = size( kspace);
n_points = prod( sz_kspace(1:(end-1)));
if sz_kspace( end) ~= 3
    error( 'kspace must be an array of 3-vectors')
end

% flattening and preparing kspace
kspace = reshape( kspace, n_points, 3);
kspace = kspace * R;
argx = .5 * kmax * kspace(:, 1) * L(1);
argy = .5 * kmax * kspace(:, 2) * L(2);
argz = .5 * kmax * kspace(:, 3) * L(3);

% calculation
kdata = ( gain*L(1)*gauss( argx) ) .* ( L(2)*gauss( argy) ) .* ( L(3)*gauss( argz) ) ./ ( exp( -i*2*pi*( .5 * kmax * kspace * D')) ); % factor of .5 is for the GE FOV

% return shaped kdata
kdata = reshape( kdata, [sz_kspace(1:(end-1)), 1]);

return
