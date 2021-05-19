% convert to forSSN_3d*.mat

%clear all;

startDir = pwd;
% [fileName, pathName, filterIndex] = uigetfile( ...
%     {'P*.mat','MAT-files (P*.mat)'}, ...
%     'Pick a file containing an ngfnRecon object...');
% 
% fprintf( 'Loading P*.mat file:  %s\n', fileName);
% fprintf( 'from directory:  %s\n', pathName);
% 
% cd( pathName);
% load( fileName);

% Arguments:
%   kdata   [np nv ns nc]   complex
%   knav    [np nvn ns nc]  complex
%   ktraj   [3 np nv ns]
%   dcf     [np nv ns]

kdata = testobj.fid;
kdata_size = size( kdata)

knav = testobj.fid_nav;
knav_size = size( knav)

ktraj = shiftdim( testobj.kspaceRad, 3);
ktraj_size = size( ktraj)

dcf = testobj.dcf;
dcf_size = size( dcf)

ktraj_nav = testobj.kspaceRad_nav;
ktraj_nav_size = size( ktraj_nav);
vars = {'kdata' 'knav' 'ktraj' 'dcf', 'ktraj_nav'};
%clear testobj;
uisave( vars);
cd( startDir);