classdef ngfnRecon
    %Matlab class for processing and reconstructing ngfnRecon data
    % outobj = ngfnRecon( varargin)
    % outobj = ngfnRecon( fidPath, viewPath, FilterIndex, SIZE, DEBUG, WARNING, ERROR, DYNAMIC and DMF)
    %
    % version for ISMRM 2021 Abstract #1184: Robust and Computationally Efficient Missing Point and Phase Estimation for Zero Echo Time (ZTE) Sequences
    % Copyright Champaign Imaging LLC
    % This work is licensed under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
    % https://creativecommons.org/licenses/by-nc-nd/4.0/
    
    %% Constant properties
    properties (Constant = true)
        version = '210414.01';  % information about the class version, constants and methods
        version_comment = 'ISMRM 2021 Abstract #1184';
    end
    
    %% Variable properties
    properties
        filterindex % serves as a datatype flag (GE, Varian, etc.)
        FLAGS       % for controlling DEBUG, WARNING, ERROR levels, DYNAMIC and DMF added
        flags_org   % flags of original reconstruction, when loading from forSSN_3d file
        
        startPath   % initial path when object constructed
        
        fidPath		% the full path to the vnmrj ".fid" directory usedngfnRecon( 'nada', 'nada', 0) for obj construction
        fidName		% name of the fid file, defaults to "fid", Pfile or silent_grad also...
        fid         % the complex fid data as fid( readout, views, channels) ...check this, now in dynamic segments CAC 200815 ***
        fid_nav     % data for navigator segments
        fid_s       % data for waspi segment
        fid_pad     % data for padding views
        
        pars		% the vnmrj parameters of interest (make an object dependent on procpar?)
        fp_procpar  % file pointer to the procpar file, could use for the Pfile also?
        
        views		% the views acquired, (x, y, z), points in gradient AD units
        viewPath
        kspaceRad	% normalized (to approximately 1 at edge) radial k-space coordinates for gridding
        kspaceRadLambda	% lambda (dipole strength) values for each coordinate
        kspaceRad_nav   % k=space coordinates of periodic navigators
        
        kspCart     % cartesian kspace values
        kspAPC      % cartesian kspace values of (dynamic) apodization correction
        
        dcf         % the density correction function (weights) corresponding to kspaceRad
        param       % dcf and gridding pipeline related parameters
        switches    % pipeline related options
        
        rftable     % SWIFT and Varian specific parameters (can these be in child object?)
        fidMainHdr
        fidBlockHdr % put these in a structure "djswift" ?
        
        SToRM       % SToRM related output
        
        % place in computed parameters structure?
        nrcvrs		% the number of channels
        
        saveDirName % name of directory to place reconstructions, logs, etc.
        savePath    % full path to place reconstructions, logs, etc.
        
        hdr         % hdr information (nifti) for the reconstructed images
        img         % complex images from all channels and frames (3d, 4d, 5d)
        hdr_sos     % hdr information (nifti) for the reconstructed rsos images
        img_sos     % rsos images for each frame (3d, 4d)
        hdr_csm     % hdr for coil sensistivity maps 
        img_csm     % coil sensistivity maps (from simulation)
        
        comment     % comment from command line, make an agument *** CAC 201207
        debug       % put debugging variables, structures, calsses here
    end
    
    %% Methods
    methods
        % constructor
        function obj = ngfnRecon( varargin)
            
            
            % Default values of FLAGS
            obj.FLAGS.SIZE = 0;     % default is autosize based on length of fid
            
            %obj.FLAGS.DEBUG = 0;   % add a command line option?
            %obj.FLAGS.DEBUG = 1;   % ngfnRecon level output/timing, log file
            obj.FLAGS.DEBUG = 2;   % method call timing, debug info in obj
            %obj.FLAGS.DEBUG = 3;   % methods internal debug, some intermediate data plots displayed and saved
            %obj.FLAGS.DEBUG = 4;   % save the obj, more plots
            %obj.FLAGS.DEBUG = 5;	% lots of stuff
            
            obj.debug.debug_level = obj.FLAGS.DEBUG;
            
            if ( obj.FLAGS.DEBUG >=3 )
                fprintf( '\n%s DEBUG ====================\n', mfilename( 'fullpath'));
            end
            
            %obj.FLAGS.WARNING = 0; % warnings off
            obj.FLAGS.WARNING = 1; % warnings on
            
            obj.FLAGS.ERROR = 1;
            
            obj.FLAGS.DYNAMIC = 0;
            
            obj.FLAGS.DMF = 0;
            
            %obj.FLAGS.SVF = 0;  % do not save any .mat data
            obj.FLAGS.SVF = 1;  % save forSSN_3d_xxxxx.mat data
            
            % load needed packages
            %pkg load signal; (octave ok
            
            % paths
            obj.startPath = getenv('NGFN_startPath');
            cd( obj.startPath);
            
            % prepare input
            if nargin < 1
                fidPath = '';
                dataDir = 'data';
                cd( dataDir);
            else
                fidPath = string( varargin{1});
                [fidPath, fidName, fidExt] = fileparts( fidPath);
                addpath( fidPath); % put on matlab path to use 'which'
                fidName = char( strcat( fidName, fidExt));
                fidPath = which ( fidName);
                [fidPath, fidName, fidExt] = fileparts( fidPath);
                obj.fidPath = char( fidPath);   obj.fidName = char( strcat( fidName, fidExt));
            end
            if nargin < 2
                viewPath = '';
                obj.pars.vieworder = 'nada';
                obj.viewPath = '';
            else
                viewPath = char( varargin{2});
                if strcmp( viewPath, 'nada')
                    obj.pars.vieworder = 'nada';
                    obj.viewPath = '';
                else
                    viewPath = which( viewPath);
                    [viewPath, viewName, EXT] = fileparts( viewPath);
                    obj.viewPath = char( viewPath); obj.pars.vieworder = char( viewName);
                end
            end
            if nargin < 3
                filterindex = 0;
            else
                filterindex = varargin{3};
            end
            
            if filterindex == 0 && strcmp( fidPath, '');
                % open a file dialog
                [obj.fidName, obj.fidPath, filterindex] = uigetfile( ...
                    {'forSSN_3d_*.mat',  'forSSN_3d*.mat data'}, ...
                    'Please select a data file to view and reconstruct');
            else
                error( 'file arguments currently not accepted, *** CAC 210514');  
            end
            
            fidName = obj.fidName;
            
            % Enable pass filterindex to methods
            obj.filterindex = filterindex;
            
            for i_arg = 4:(nargin)
                field = {'SIZE', 'DEBUG', 'WARNING', 'ERROR', 'DYNAMIC', 'DMF'};
                obj.FLAGS = setfield( obj.FLAGS, field{i_arg-3}, varargin{i_arg});
            end
            
            if ( obj.FLAGS.DEBUG >= 1 ); tstart = tic; end % start timer
            if ( obj.FLAGS.WARNING >= 1 ); warning('ON'); else warning('OFF'); end % set warning
            
            % create the saved image directory, add some error checking?
            switch  filterindex
                case 0  % cancelled
                    fprintf( '\n'); warning( 'Cancelled by user');
                    % cleanup
                    cd( obj.startPath);
                    % total elapsed time
                    if ( obj.FLAGS.DEBUG >= 1 )
                        fprintf( 'ngfnRecon Total ');    toc( tstart)
                    end
                    % close diary
                    %if ~( filterindex == 1 || filterindex == 11 )
                    if ~( filterindex == 1 )
                        if ( obj.FLAGS.DEBUG >= 1 )
                            diary off;
                        end
                    end
                    return;
                case 1 % forSSN_3d mat file
                    saveDirName = '_recon_';
                otherwise
                    fprintf( '\n'); error( 'unknown case for filterindex, *** CAC 180227');
            end
            
            reply = input( 'How many missing points to estimate? [0]:');
            if isempty( reply)
                Missing_Points = 0;
            else
                Missing_Points = reply;
            end
            
            STR = ''; STR = input('Enter a comment or return for none: ','s');
            if ~isempty( STR); saveDirName = strcat( saveDirName, STR, '_'); end;
            
            
            
            % check if exists
            if exist( obj.fidPath, 'dir')
                %everything ok
            else
                % code for ERROR_FLAG ***
                fprintf( '\n'); error( 'directory does not exist');
            end
            
            % Create working directory
            %   savedir paths, currently create in same directory as data

                [PATHSTR, fidName, EXT] = fileparts( obj.fidName);  % How to better pick the name?????????, ********* CAC 190524 *********
                
                obj.saveDirName = strcat( fidName, saveDirName, datestr( now, 30));
                obj.savePath = strcat( obj.fidPath, '/', obj.saveDirName);
                mkdir( obj.fidPath, obj.saveDirName);
                cd( obj.savePath);
                if ( obj.FLAGS.DEBUG >= 1 )
                    diary_file = strcat( obj.saveDirName, '.log');
                    diary( diary_file);
                end

            
            if ( obj.FLAGS.DEBUG >= 3 ); fprintf( '\n%s DEBUG obj.fidPath: %s', mfilename, obj.fidPath); end
            if ( obj.FLAGS.DEBUG >= 3 ); fprintf( '\n%s DEBUG obj.fidName: %s', mfilename, obj.fidName); end
            if ( obj.FLAGS.DEBUG >= 3 ); fprintf( '\n%s DEBUG obj.savePath: %s\n\n', mfilename, obj.savePath); end
            
            if ( obj.FLAGS.DEBUG >= 1 )
                objver = obj.version;  dtstr = datestr( now);
                fprintf( '===== ngfnRecon version = %s =====\n', objver);
                fprintf( 'version_comment: %s\n', obj.version_comment);
                fprintf( 'runtime: %s\n', dtstr);
                fprintf( 'FLAGS: SIZE = %d, ', obj.FLAGS.SIZE);
                fprintf( 'DEBUG = %d, ', obj.FLAGS.DEBUG);
                fprintf( 'WARNING = %d, ', obj.FLAGS.WARNING);
                fprintf( 'ERROR = %d, ', obj.FLAGS.ERROR);
                fprintf( 'DYNAMIC = %d, ', obj.FLAGS.DYNAMIC);
                fprintf( 'DMF = %d\n', obj.FLAGS.DMF);
                % make a string array for action description, CAC 181218 ***
                fprintf( 'carrying out action %d ...\n', obj.filterindex);
            end
            
            if ( 3 > obj.FLAGS.DEBUG >= 1); fprintf( 'MATLAB Version: %s\n', version); end % print matlab version
            if ( obj.FLAGS.DEBUG >= 3 ); ver; end % print detailed matlab and system info
            
            % show license and user info
            [status, result] = unix( 'whoami'); fprintf( 'Username: %s', result);
            [status, result] = unix( 'hostname'); fprintf( 'Hostname: %s', result);
            
            maxNumCompThreads( 'automatic'); nThreadsAuto = maxNumCompThreads;
            maxNumCompThreads( 2*nThreadsAuto); nThreads = maxNumCompThreads;
            if ( obj.FLAGS.DEBUG >= 1 ); fprintf( 'nThreads = %d\n', nThreads); end % print number of computation threads
            
            [gpu_status, gpu_result] = system( 'nvidia-smi -L'); % check for gpus, Driver Version
            if ( obj.FLAGS.DEBUG >= 1 ) & ~gpu_status; fprintf( 'GPU Information...\n%s', gpu_result);  end %gpu_result(1:(end-1)) removes newline or return
            cd( obj.savePath);
            
            cd( obj.startPath);
            [git_status, result] = system( 'git status'); % check for git directory
            if ( obj.FLAGS.DEBUG >= 1 ) & ~git_status; fprintf( 'Git Information: '); system( 'git describe --all --dirty --long'); end
            if ( obj.FLAGS.DEBUG >= 2 ) & ~git_status; system( 'git status --porcelain'); end
            if ( obj.FLAGS.DEBUG >= 4 ) & ~git_status; system( 'git status -vv'); end
            
            
            % Menu actions
            %   put into separate calls/files, *** CAC 190223
            
            % forSSN_3D_*.mat file
            if filterindex == 1
                if strncmp( obj.fidName, 'forSSN_3d_', 10) %%% eventually make this a case sctucture, *** CAC 180110
                    obj = obj.get_forSSN_3d;
                    obj.pars.rcvr_gate = Missing_Points;
                    obj = obj.doAutophase;
                    obj = obj.doKdata;
                    obj = obj.doImages;
                end
            end
            
            % total elapsed time
            if ( obj.FLAGS.DEBUG >= 1 )
                fprintf( 'ngfnRecon Total ');    toc( tstart)
            end
            
            %             if ( obj.FLAGS.DEBUG >= 2 )
            %                 obj.comment = input( 'Comment for this reconstruction: ', 's');
            %             end
            
            if filterindex ~= 2
                if ( obj.FLAGS.DEBUG >= 1 )
                    diary off;
                end
            end
            
            cd( obj.startPath);
            
        end
        
    end
    
    %% External methods
    methods     % in external files
        obj = get_forSSN_3d( obj);       % loads a forSSN_3d mat file
        
        obj = doDCF( obj);               % calculates the dcf
        obj = doKdata( obj);             % grids radial k-space and corrleated fid values to cartesian k-space values
        obj = doImages( obj);            % makes images and image series from kspCart and kspAPC
        obj = doAutophase( obj);         % phase and missing points estimation
    end
    
end
