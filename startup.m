% startup.m

% for docker container and git cloning, set up path on first run
ngfnPathDef = '../ngfn_simulation_recon/pathdef.m';
ngfnPath = genpath( '../ngfn_simulation_recon');
if exist( ngfnPathDef) == 0
    fprintf( '\n'); warning( 'Creating initial %s\n', ngfnPathDef);
    addpath( pathdef); addpath( ngfnPath); savepath( ngfnPathDef);
end

% show version and packages
ver; fprintf( '-----------------------------------------------------------------------------------------------------\n');

% cpu info
cpuInfo = cpuinfo;
fprintf( 'CPU: %s\n', cpuInfo.Name)
fprintf( 'Clock: %s  Cores: %d  Cache: %s\n', cpuInfo.Clock, cpuInfo.NumProcessors, cpuInfo.Cache);
fprintf( 'OSType: %s  OSVersion: %s\n', cpuInfo.OSType, cpuInfo.OSVersion);
[status, result] = unix( 'cat /proc/meminfo | grep MemFree:'); fprintf( '%s', result);
fprintf( '-----------------------------------------------------------------------------------------------------\n');

% gpu info
[Ngpu] = gpuDeviceCount; fprintf( 'GPU device count: %d\n', Ngpu);
for idx = 1:Ngpu;
    gpuInfo( idx) = gpuDevice( idx);
    fprintf( 'GPU device: %d  Name: %s  MultiprocessorCount: %d  AvailableMemory: %d MB\n', ...
        gpuInfo( idx).Index, gpuInfo( idx).Name, gpuInfo( idx).MultiprocessorCount, floor( gpuInfo( idx).AvailableMemory/1024/1024));
end
fprintf( '-----------------------------------------------------------------------------------------------------\n');

% show license and user info
[status, result] = unix( 'hostname'); fprintf( 'Hostname: %s', result);
[status, result] = unix( 'whoami'); fprintf( 'Username: %s', result);

% python version
pythonEnv = pyenv; pyVer = pythonEnv.Version; fprintf( 'Python: v%s\n', pyVer);

% start time
fprintf( 'Start Time: %s\n', datestr(now));
fprintf( '-----------------------------------------------------------------------------------------------------\n');
fprintf( 'to get started type:  testobj = ngfnRecon;\n');  
fprintf( '-----------------------------------------------------------------------------------------------------\n');