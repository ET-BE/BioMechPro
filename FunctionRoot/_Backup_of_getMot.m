function [] = getMot(C3Ddata,filename,varargin)
%% getMot
% Create .mot file with external loads (ground reaction forces and moments)
% for use of point force application in OpenSim.
% 
% INPUT)
% C3Ddata : Data structure containing at least the fields:
% [..TODO..]
% 
% filename : Path and name of the .mot file to be created.
% If no path is provided the file is saved in the current folder.
% 
% permvec : permutation vector to put data in OpenSim coordinates
% (x: walking direction, z: to the right, y: upward)
% Default is [1 2 3] (no permutation)
% 
% extraXLD : string or cell of strings, containing label names in
% OtherDataLabel. Corresponding columns of OtherData will also be written
% to the file as one-dimensional external loads (see extraXLDdim).
% 
% extraXLDP : 3 element vector, or Nx3 matrix, specifying point data that
% indicates where each extraXLD is applied. For the ground reaction force
% this is equivalent to the COP (without torque) or the origin of the force
% plate (with torque). Input extraXLDP must have 1 row for every extraXLD.
% It will default to 0 if it is not specified.
% Example: Suppose you want to apply an additional external load to the COM 
% of the pelvis. The origin of the pelvis segment in OpenSim is not its COM.
% Consequently, you need to apply a point force to the pelvis, with the 
% point expressed in pelvis, specifying the distance (x,y,z) of the COM 
% relative to the pelvis origin. The force itself might be expressed in any 
% other coordinate system (e.g. ground).
% 
% extraXLDdim : numerical vector, containing a dimension index 1,2 or 3
% specifying in which OpenSim dimension the extraXLD channels are stored.
% The numbers 1-3 give suffix 'x','y', and 'z' respectively.
% The two other not specified dimensions will contain only zeros.
% The number of elements in extraXLDdim must be the same as the number of 
% channels specified by extraXLD. If extraXLDdim is not specified, it will
% default to 1 for every channel in extraXLD.
% 
% 
% OUTPUT)
% No direct output
% A .mot file is created in the destination provided by filename
% 
% NOTES)
% This file does not do any processing of the force data. 
% It ports force data to a .mot file, where the number of samples written
% corresponds with the number of marker samples for which force data is
% available. It is assumed that the force data has AT LEAST the sample
% frequency of that of the marker data. If the force data has a higher
% sample frequency, some samples will be discarded.
% 
% A dual plate is assumed, one for each foot
% Channel order is assumed:
% Forces Left (1-3)
% Moments Left (4-6)
% Forces Right (7-9)
% Moments Right (10-12)
% 
% In OpenSim, you should apply a POINT FORCE to the foot, and add an
% additional TORQUE. Both the point force and the torque should be 
% expressed in GROUND coordinates. 
% If MoCap and force plate have the same origin, then these point coordinates are
% [0 0 0], as both force and moment originate from the origin of the force plate.
% In other words, a force vector and a moment originating from the plate's 
% origin is the same as a force vector originating from a point outside 
% the plate's origin (i.e. in the COP).
% 
% All external forces MUST have 3 components with a common prefix in the
% header, and end the name on either x,y or z. 
% This means you cannot specify individual channels as forces or moments 
% applied to the model. So for every channel specified in extraXLD, two 
% additional zero channels with the same prefix will be created.

% Mark Vlutters - September 2015 - Enschede

%% Settings

% Set vertical force threshold
% fthresh = 20; % default

%% Check input

if nargin == 2
    permvec = [1 2 3];
    extraXLD = {};
    extraXLDdim = [];
    extraXLDP = [];
elseif nargin == 3
    permvec = varargin{1};
    if (size(permvec,1) == 3) && (size(permvec,2) == 1)
        permvec = permvec';
    elseif sum(size(permvec) - [1 3]) ~= 0
        error('getMot:permvec','Input parameter permvec must be a 1x3 vector.');
    end
    extraXLD = {};
    extraXLDdim = [];
    extraXLDP = [];
elseif nargin == 4
    permvec = varargin{1};
    if (size(permvec,1) == 3) && (size(permvec,2) == 1)
        permvec = permvec';
    elseif sum(size(permvec) - [1 3]) ~= 0
        error('getMot:permvec','Input parameter permvec must be a 1x3 vector.');
    end
    extraXLD = varargin{2};
    if ischar(extraXLD) 
        extraXLD = {extraXLD}; % String to cell
    end
    extraXLDdim = ones(1,length(extraXLD));
    extraXLDP = zeros(length(extraXLD),3);
elseif nargin == 5
    permvec = varargin{1};
    if (size(permvec,1) == 3) && (size(permvec,2) == 1)
        permvec = permvec';
    elseif sum(size(permvec) - [1 3]) ~= 0
        error('getMot:permvec','Input parameter permvec must be a 1x3 vector.');
    end
    
    extraXLD = varargin{2};
    if ischar(extraXLD) 
        extraXLD = {extraXLD}; % String to cell
    end
    
    extraXLDdim = varargin{3};
    if numel(extraXLDdim) ~= numel(extraXLD)
        error('getMot:extraXLDdimEl','Number of elements in extraXLD must be the same as in extraXLDdim');
    elseif any((extraXLDdim<1)|(extraXLDdim>3))
        error('getMot:extraXLDdim','extraXLDdim may not contain elements outside the range 1:3.');
    elseif size(extraXLDdim,2) == 1
        extraXLDdim = extraXLDdim';
    end
    
elseif nargin > 5
    error('getMot:inputs','Too many input arguments');
end

% Check trialname
if ~ischar(filename)
    error('getMot:trialname','Input trialname must be a string');
end

%% Collect some info

% If full path is supplied, take last part for name inside mot file
if ~isempty(strfind(filename,'\'))
    foo = strfind(filename,'\');
    infilename = filename(foo(end)+1:end);
    pathname = filename(1:foo(end));
else
    infilename = filename;
    pathname = pwd;
end

% Check if the extra fields exist in OtherData
if isfield(C3Ddata,'Other') && ~isempty(extraXLD)
    for ixld = 1:length(extraXLD)
        addIdx(ixld) = find( strcmpi(C3Ddata.Other.OtherDataLabel,extraXLD{ixld}) ); % NOTE: Assumed 2D data, and there is no check for non-existing name
    end
elseif ~isempty(extraXLD)
    warning('getMot:extraXLD','No Other field in structure. Cannot add additional external loads.');
end

% Check if sync indices exist for force data
if ~isfield(C3Ddata.Marker,'MarkerSyncIdx') || ~isfield(C3Ddata.Force,'ForceSyncIdx')
    warning('getMot:syncIdx',['No sync indices found in' infilename '. Assuming synchronized marker-force data']);
    
    markerSyncIdx = [1 size(C3Ddata.Marker.MarkerData,1)];
    forceSyncIdx = [1 size(C3Ddata.Force.ForceData,1)];
else
    markerSyncIdx = C3Ddata.Marker.MarkerSyncIdx;
    forceSyncIdx = C3Ddata.Force.ForceSyncIdx;
end

% Check if sync indices exist for additional data
if ~isempty(extraXLD)
    if ~isfield(C3Ddata.Other,'OtherSyncIdx')
        otherSyncIdx = [1 size(C3Ddata.Other.OtherData,1)];
    else
        otherSyncIdx = C3Ddata.Other.OtherSyncIdx;
    end
end

% Frame rates
markerFrameRate = C3Ddata.Marker.MarkerFrameRate;
% forceFrameRate = C3Ddata.Force.ForceFrameRate;

% Data column header
columnNames = {'time' ...
    'l_ground_force_vx' 'l_ground_force_vy' 'l_ground_force_vz' ... % Left GRF Vector in specific body CRF
    'l_ground_torque_x' 'l_ground_torque_y' 'l_ground_torque_z' ... % Left Moments
    'ground_force_vx' 'ground_force_vy' 'ground_force_vz' ... % Right GRF Vector in specific body CRF
    'ground_torque_x' 'ground_torque_y' 'ground_torque_z' ... % Right Moments
    'l_ground_force_px' 'l_ground_force_py' 'l_ground_force_pz' ... % Left COP
    'ground_force_px' 'ground_force_py' 'ground_force_pz' }; % Right COP

% Add additional external loads, if any
if ~isempty(extraXLD)
    for ixld = 1:length(extraXLD);
        columnNames(end+1:end+3) = {[extraXLD{ixld} '_x'] , [extraXLD{ixld} '_y'] , [extraXLD{ixld} '_z']};
    end
    
    if ~isempty(extraXLDPt)
        columnNames(end+1:end+3) = {[extraXLD{ixld} '_px'] , [extraXLD{ixld} '_py'] , [extraXLD{ixld} '_pz']};
    end
end

% Data size
nRows = 1 + markerSyncIdx(end) - markerSyncIdx(1);
nColumns = length(columnNames);

%% Get number of force samples equal to number of marker samples

nMarkSmpl = markerSyncIdx(end) - markerSyncIdx(1) + 1;
forceIdx = round(linspace(forceSyncIdx(1),forceSyncIdx(end),nMarkSmpl));

if ~isempty(extraXLD)
    otherIdx = round(linspace(otherSyncIdx(1),otherSyncIdx(end),nMarkSmpl));
end

%% Collect FORCE, MOMENT and COP data to write
% (currently, no COP data is put into the .mot file)

% Collect force and moment data
permall = [permvec permvec+3 permvec+6 permvec+9];  % Assumed 12 channel MGRF
forceData = C3Ddata.Force.ForceData(forceIdx,permall);

if ~isempty(extraXLD)
    % NOTE: permvec is not applied here!
    otherData = zeros(length(otherIdx),3*length(addIdx));
    otherData(:,extraXLDdim + (0:3:3*length(addIdx)-1) ) = C3Ddata.Other.OtherData(otherIdx,addIdx);
else
    otherData = [];
end

% COP data (zeros, as we apply a force AND a torque, 
% unless the forceplate CRF has an offset wrt the MoCap CRF)
% TODO: allow for time varying plate offsets (timeseries in structure)
if ~isfield(C3Ddata.Force,'ForcePlateOffset')
    copData = zeros(size(forceData,1),6);
else
    copData = repmat(C3Ddata.Force.ForcePlateOffset,[size(forceData,1) 1]);
    copData = copData(:,[permvec permvec+3]);
end

% Add time column
writeData = [(0:size(forceData,1)-1)'./markerFrameRate forceData copData otherData];

%% Create file and write header

fid = fopen([filename 'XLD.mot'],'w');
if fid == -1
    error('getMot:FileID',['Cannot open ' infilename '.mot for writing. It might be in use.']);
end

% General header (SIMM header only)
fprintf(fid, [infilename '\n' ...
    'version=1\n' ...
	'nRows=' num2str(nRows) '\n' ...
    'nColumns=' num2str(nColumns) '\n' ...
    'inDegrees=no\n' ...
    'endheader\n']);

% Column headers
for iCol = 1:nColumns
    fprintf(fid, [columnNames{iCol} '\t']);
end
fprintf(fid, '\n');

%% Write data

writeStr = regexprep(mat2str(writeData),{'[',']',' ',';'},{'','','\t','\n'});
fprintf(fid,writeStr);

%% Clean up

fclose(fid);

disp([infilename 'XLD.mot created in ' pathname]);

end