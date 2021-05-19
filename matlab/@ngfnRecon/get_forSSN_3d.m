function obj = get_forSSN_3d(obj)
%ngfnRecon.get_forSSN_3d
% load data from forSSN_3d*.mat file
% 
% C. Corum 5/13/2021
% Copyright Champaign Imaging LLC
% This work is licensed under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
% https://creativecommons.org/licenses/by-nc-nd/4.0/


DEBUG_FLAG = obj.FLAGS.DEBUG;
%DEBUG_FLAG = 3;

if ( DEBUG_FLAG >=3 )
    fprintf( '\n%s DEBUG ====================\n', mfilename( 'fullpath'));
end

if ( DEBUG_FLAG >= 2 )
    gfssn_timer = tic;    fprintf('Loading .mat data from %s, ', obj.fidName);
end

directory       = obj.fidPath;
fidPath         = strcat( directory, '/', obj.fidName);

cd( directory)
load( fidPath)

obj.fid = kdata;

obj.fid_nav = knav;

obj.kspaceRad = squeeze( shiftdim( ktraj, 1)); 
                
obj.dcf = dcf;
                
obj.kspaceRad_nav = ktraj_nav; % check this
                
obj.img_csm = csm;
                
obj.pars = pars;

obj.views = views;            

obj.switches = switches;

obj.param = param;

obj.flags_org = flags_org;

%return % **********************  WORKING HERE ***************************** CAC 210513

% write switches,
obj.pars.rcvrs_recon = obj.pars.rcvrs;     % all channels

% recievers subset
%obj.pars.rcvrs_recon  = obj.pars.rcvrs * 0;    % DO THIS FIRST!
%obj.pars.rcvrs_recon(1) = 1;                   % first channel only
%obj.pars.rcvrs_recon(35) = 1;                  % channel 35 only
%ch = [34:48]; obj.pars.rcvrs_recon(ch) = 1; % channels in list only (use numbers for virtual coil groups? *** CAC 210329)

obj.pars.djNch_recon = sum( obj.pars.rcvrs_recon);

% sizes and fix for channels = 1, volumes = 1;
size_fid = size( obj.fid);
if ndims( obj.fid) < 3
    size_fid(3) = 1;
end
views_fid = size_fid(2) * size_fid(3);
ndims_fid = ndims( obj.fid);
if ndims_fid < 4
    size_fid(4) = 1;
end

% sizes and fix for channels = 1, volumes = 1;
size_fid_nav = size( obj.fid_nav);
if ndims( obj.fid_nav) < 3
    size_fid_nav(3) = 1;
end
views_fid_nav = size_fid_nav(2) * size_fid_nav(3);
ndims_fid_nav = ndims( obj.fid_nav);
if ndims_fid_nav < 4
    size_fid_nav(4) = 1;
end


if ( DEBUG_FLAG >= 2 )
    fprintf('Loaded %d highres and %d navigator views from %d channel(s), ', views_fid, views_fid_nav, size_fid(4));  toc( gfssn_timer);
end

return
