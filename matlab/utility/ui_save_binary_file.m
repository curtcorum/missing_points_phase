
% ui_save_binary_file
%   Curt Corum, 210422

FilterSpec = '*.csg';
Title = 'Save to a binary file';
[FileName, PathName] = uiputfile( FilterSpec, Title);
FilePath = fullfile( PathName, FileName);
%fprintf( '%s\r', FilePath);

fidtosave = testobj.fid(:);

write_complex_binary (fidtosave, FilePath);

% fileID = fopen( FilePath, 'w');
% count = fwrite( fileID, fidtosave, 'float32')
% st = fclose( fileID)
% plot(RawData);


