function pArray = pointArray( varargin)
%function pArray = pointArray( N, S, R)
% creates a 3d array of 3-vectors
%
% defaults: N = [1 1 1]; S = [1 1 1], R = eye(3);
%
% Curt Corum, Champaign Imaging LLC, 5/16/2019

% defaults:
N = [3 3 3]; S = [1 1 1]; R = eye(3);

% varargin handling
if nargin > 0; N = varargin{1}; end
if nargin > 1; S = varargin{2}; end
if nargin > 2; R = varargin{3}; end
if nargin > 3; error('too many arguments'); end

% input checking

% prepare variables
Nx = N(1); Ny = N(2); Nz = N(3);
Cx = (Nx+1)/2; Cy = (Ny+1)/2; Cz = (Nz+1)/2;
Sx = S(1); Sy = S(2); Sz = S(3);

% fill in vectors for array
pArray = zeros( Nx, Ny, Nz, 3);
parfor idx = 1:Nx
    for idy = 1:Ny
        for idz = 1:Nz
            pArray(idx, idy, idz, :) = [(idx-Cx)*Sx (idy-Cy)*Sy (idz-Cz)*Sz]*R';
        end
    end
end

return
