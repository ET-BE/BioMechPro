function [] = getTrc(C3Ddata,filename,varargin)
%% getTrc
% Get .trc file from data, for example for use in OpenSim
% This file is required for scaling a general simbody model to subject characteristics.
% 
% For more information on .trc files see:
% 
% http://simtk-confluence.stanford.edu:8080/display/OpenSim/Marker+%28.trc%29+Files
% 
% For more information on inverse kinematics and dynamics see:
% 
% http://simtk-confluence.stanford.edu:8080/display/OpenSim/Tutorial+3+-+Scaling,+Inverse+Kinematics,+and+Inverse+Dynamics
% 
% INPUT)
% C3Ddata : Data structure containing at least the fields:
% [...]
% 
% trialname : Path and name of the .trc file to be created.
% If no path is provided the file is saved in the current folder.
% 
% dataflag : flag to specify which data to write to the .trc file
% 1: marker data only
% 2: probe data only
% 3: marker and probe data
% Default value is 1
% 
% permvec : permutation vector to put data in OpenSim coordinates
% (x: walking direction, z: to the right, y: upward)
% Default is [1 2 3] (no permutation)
% 
% OUTPUT)
% No direct output
% A .trc file is created in the destination provided by trialname

% Mark Vlutters - July 2015 - Enschede

%% Check input

if nargin == 2
    dataflag = 1;
    permvec = [1 2 3];
elseif nargin == 3
    dataflag = varargin{1};
    permvec = [1 2 3];
elseif nargin == 4
    dataflag = varargin{1};
    permvec = varargin{2};
elseif nargin > 4
    error('getTrc:inputs','Too many input arguments');
end

% Check if fields are compatible with the flag
if ( (dataflag == 1) && ~isfield(C3Ddata.Marker,'MarkerData') ) || ...
        ( (dataflag == 2) && ~isfield(C3Ddata.Marker,'ProbedData')) || ...
        ( (dataflag == 3) && (~isfield(C3Ddata.Marker,'MarkerData') || ~isfield(C3Ddata.Marker,'ProbedData')) )
    error('getTrc:datafield','MarkerData or ProbedData field unavailable in C3Ddata structure');
end

% Check trialname
if ~ischar(filename)
    error('getTrc:trialname','Input trialname must be a string');
end

%% Collect some info
markerFrameRate = C3Ddata.Marker.MarkerFrameRate;

if dataflag == 1
    nMarker = size(C3Ddata.Marker.MarkerData,2);
elseif dataflag == 2
    nMarker = size(C3Ddata.Marker.ProbedData,2);
elseif (dataflag == 3)
    nMarker = size(C3Ddata.Marker.MarkerData,2) + size(C3Ddata.Marker.ProbedData,2);
end

% If full path is supplied, take last part for name inside trc file
if ~isempty(strfind(filename,'\'))
    foo = strfind(filename,'\');
    infilename = filename(foo(end)+1:end);
    pathname = filename(1:foo(end));
else
    infilename = filename;
    pathname = pwd;
end

% Check if sync indices exist
if ~isfield(C3Ddata.Marker,'MarkerSyncIdx') || ~isfield(C3Ddata.Force,'ForceSyncIdx')
    warning('getTrc:syncIdx',['No sync indices found in ' infilename '. Assuming synchronized marker-force data (if any)']);
    
    markerSyncIdx = [1 size(C3Ddata.Marker.MarkerData,1)];
else
    markerSyncIdx = C3Ddata.Marker.MarkerSyncIdx;
end

nFrames = 1 + markerSyncIdx(end) - markerSyncIdx(1);

%% Collect header to write

if dataflag == 1
    writeHeader = C3Ddata.Marker.MarkerDataLabel;
elseif dataflag == 2
    writeHeader = C3Ddata.Marker.ProbedDataLabel;
elseif dataflag == 3
    % Note: write probe data first, probably used most
    writeHeader = [C3Ddata.Marker.ProbedDataLabel C3Ddata.Marker.MarkerDataLabel];
end

%% Collect data to write

if dataflag == 1

    markerData = C3Ddata.Marker.MarkerData(markerSyncIdx(1):markerSyncIdx(end),:,permvec);
    
    msiz = size(markerData);
    
    writeData = [...
        (1:msiz(1))'  , ...
        (0:msiz(1)-1)'./markerFrameRate , ...
        reshape(permute(markerData,[1 3 2]),[msiz(1) msiz(2).*msiz(3)]) ];
    
elseif dataflag == 2
    
    probedData = C3Ddata.Marker.ProbedData(markerSyncIdx(1):markerSyncIdx(end),:,permvec);
    
    psiz = size(probedData);
    
    writeData = [...
        (1:psiz(1))' , ...
        (0:psiz(1)-1)'./markerFrameRate , ...
        reshape(permute(probedData,[1 3 2]),[psiz(1) psiz(2).*psiz(3)]) ];
    
elseif dataflag == 3
    
    markerData = C3Ddata.Marker.MarkerData(markerSyncIdx(1):markerSyncIdx(end),:,permvec);
    probedData = C3Ddata.Marker.ProbedData(markerSyncIdx(1):markerSyncIdx(end),:,permvec);
    
    msiz = size(markerData);
    psiz = size(probedData);

    % Note: write probe data first, probably used most
    writeData = [...
        (1:msiz(1))' , ...
        (0:msiz(1)-1)'./markerFrameRate , ...
        reshape(permute(probedData,[1 3 2]),[psiz(1) psiz(2).*psiz(3)]) , ...
        reshape(permute(markerData,[1 3 2]),[msiz(1) msiz(2).*msiz(3)]) ];
    
end

%% Create file and headers
% File will contain all markers and probe positions
% Not all might have to be associated with the scaling of the model

% File
fid = fopen([filename '.trc'],'w'); % note: w also discards existing content, if any
if fid == -1
    error('getTrc:FileID',['Cannot open ' infilename '.trc for writing. It might be in use.']);
end

% General header
fprintf(fid,['PathFileType\t4\t(X/Y/Z)\t' infilename '\n' ...
    'DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n' ...
    num2str(markerFrameRate) '\t' ... % DataRate
    num2str(markerFrameRate) '\t' ... % CameraRate
    num2str(nFrames) '\t'  ... % NumFrames
    num2str(nMarker) '\t' ... % NumMarkers
    'm\t' ... % Units
    num2str(markerFrameRate) '\t' ... % OrigDataRate
    '1\t' ... % OrigDataStartFrae
    num2str(nFrames) '\n']); % OrigNumFrames

% Marker header
fprintf(fid,'Frame#\tTime\t');
for ihead = 1:length(writeHeader)
    if ihead ~= length(writeHeader)
        fprintf(fid,[writeHeader{ihead} '\t\t\t']);
    else
        fprintf(fid,[writeHeader{ihead} '\n']);
    end
end

% Dimension header
fprintf(fid,'\t\t');
for ihead = 1:length(writeHeader)
    if ihead ~= length(writeHeader)
        fprintf(fid,['X' num2str(ihead) '\tY'  num2str(ihead) '\tZ' num2str(ihead) '\t']);
    else
        fprintf(fid,['X' num2str(ihead) '\tY'  num2str(ihead) '\tZ' num2str(ihead) '\n\n']);
    end
end

%% Write data

writeStr = regexprep(mat2str(writeData),{'[',']',' ',';'},{'','','\t','\n'});
fprintf(fid,writeStr);

%% Clean up

fclose(fid);

disp([infilename '.trc created in ' pathname]);

end