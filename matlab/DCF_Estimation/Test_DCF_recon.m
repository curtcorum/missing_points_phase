% Test_DCF_recon
% loads ryan_kspace data and performs gridding reconstruction with
% iterative DCF using Rajesh's matlab and mex code
% Curt Corum, 3/18/2011

% get start time
startTime = clock;

% file dialog and change directory
startdir=pwd;
cd( '~/vnmrsys/data');
[ rkFile, rkPath] = uigetfile( '', 'Please choose a ryan kspace file');
cd( rkPath);

%rkFile='ryan_kspace'

% comment inserted here

% load kspace data and coordinates
[ kdata, kspace] = Read_Ryan_Kspace( rkFile);
cd(startdir);

% set up PARAM structure
PARAM.kernelwidth=4;
PARAM.kernelsampling=500;
PARAM.freq_scale = [1 1 1];
PARAM.oversampling = 2;
PARAM.red_matrix = 128;

% other constants
nIter=2;

% calculate DCF
[dP,d0] = Calculate_DCF_ModPipe(kspace,PARAM,nIter);

% save the DCF as ryan_kspace file
cd( rkPath);
rkFileSave = 'ryan_kspace_dcf';
Write_Ryan_Kspace( rkFileSave, dP, kspace);
cd(startdir);
fprintf( '\nryan_kspace_dcf saved\n');

% get end time
endTime = clock;

% elapsed time
elapsedTime = endTime - startTime

% exit matlab
% exit;



