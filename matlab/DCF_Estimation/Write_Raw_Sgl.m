function  Write_Raw_Sgl( magFile, phaseFile, dirName, imageComplex)
% function  Write_Raw_Sgl( magFile, phaseFile, dirName, imageComplex)
%   Writes a Mag and Phase raw files in dir, switches re and im of imageComplex
%   C. Corum, 8/21/2012
%
%   CAC 190919 - if filename 'NONE', will skip that file

rSize = size( imageComplex);

% create dir
[success, message] = mkdir( dirName);

% change to the new dir
startdir = pwd;
cd( dirName);

% open the files, big endian format, for appending data
if ~strcmp( magFile, 'NONE')
    [ magFid] = fopen( magFile, 'a+', 'ieee-be');
    imageMag = abs( imageComplex); % reformat data
    fwCount = fwrite( magFid, imageMag, 'single'); % write values
    if fwCount < rSize; error('fwCount < rSize'); end
    stDebugr = fclose( magFid); %close the file
end

if ~strcmp( phaseFile, 'NONE')
    [ phaseFid] = fopen( phaseFile, 'a+', 'ieee-be');
    imagePhase = angle( imageComplex); % reformat data
    fwCount = fwrite( phaseFid, imagePhase, 'single'); % write values
    if fwCount < rSize; warning('fwCount < rSize'); end
    stDebugi = fclose( phaseFid); %close the file
end

% other formats? CAC 190919
%imageRe = real(imageComplex);
%imageIm = imag(imageComplex);

cd( startdir);

%figure; plot( rkRe( 1:255) );
%figure; plot( rkIm( 1:255) );
%figure; plot( rkZ( 1:255) );
%plot3( rkX,  rkY, rkZ);
%plot(kdata)

end

