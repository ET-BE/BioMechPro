function [] = getMot_mod(C3Ddata,filename,varargin)
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
% You do not have to supply the .mot extension in the filename.
% 
% permvec : permutation vector to put data in OpenSim coordinates.
% (x: walking direction, z: to the right, y: upward)
% Default is [1 2 3] (no permutation)
% This vector is applied to both the ForceData in the data structure, as
% well as the extraXLD that is optionally specified.
% 
% extraXLD : string or cell of strings, containing label names in
% OtherDataLabel. Corresponding columns of OtherData will also be written
% to the file as external loads.
% 
% extraXLDdim : numerical vector, containing one or more scalar dimension 
% indices, specifying in which OpenSim dimension the extraXLD channels are 
% stored. The numbers 1:3 give suffix 'x','y', and 'z' respectively.
% Use 4:6 for the next group of variables, 7:9 for the next, etc.
% Not specified dimensions will contains zeroes.
% The number of elements in extraXLDdim must be the same as the number of 
% channels specified by extraXLD. If extraXLDdim is not specified, each
% channel specified by extraXLD will be stored in the x dimension.
% 
% OUTPUT)
% No direct output
% A .mot file is created in the destination provided by filename
% 
% NOTES)
% This function ports force data to a .mot file, where the number of 
% samples written corresponds with the number of marker samples for which 
% force data is available. 
% It is assumed that the force data has AT LEAST the sample frequency of 
% that of the marker data. If the force data has a higher sample frequency, 
% some samples will be discarded.
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
% applied to the model. So channels specified by extraXLD will also be augmented 
% to three-channel data, with empty channels containing only zeroes.
% 
% You can use extraXLD to write any kind of data in OtherData to the 
% .mot file. Not necessarily an external load, but also point data.
% 
% Example:
% >> getMot(C3Ddata,'MyFileName',[2 3 1],{'ChanA','ChanB','ChanC'},[1 2 4]);
% Creates six additional columns in the .mot file, with headers:
% ChanA_x, ChanA_y, ChanA_z, ChanC_x, ChanC_y, ChanC_z
% Containing the following data:
% ChanA_x : data from ChanB (because of permvec)
% ChanA_y : zeroes
% ChanA_z : data from ChanA
% ChanC_x : zeroes
% ChanC_x : zeroes
% ChanC_x : data from ChanC

% Mark Vlutters - September 2015 - Enschede

%% Settings

% Set vertical force threshold
fthresh = 20; % default

%% Check input

if nargin > 5
    error('getMot:inputs','Too many input arguments');
end

% Defaults
permvec = [1 2 3];
extraXLD = {};
extraXLDdim = [];

if nargin >= 3
    
    permvec = varargin{1};
    if isempty(permvec)
        permvec = [1 2 3];
    elseif (size(permvec,1) == 3) && (size(permvec,2) == 1)
        permvec = permvec';
    elseif sum(size(permvec) - [1 3]) ~= 0
        error('getMot:permvec','Input parameter permvec must be a 1x3 vector.');
    end
    
end    
if nargin >= 4
    
    extraXLD = varargin{2};
    if ischar(extraXLD) 
        extraXLD = {extraXLD}; % String to cell
    end
    
end
if nargin == 5

    extraXLDdim = varargin{3};
    dimcount = hist(extraXLDdim,max(extraXLDdim));
    if numel(extraXLDdim) ~= numel(extraXLD)
        error('getMot:extraXLDdimEl','Number of elements in extraXLD must be the same as in extraXLDdim');
    elseif any(dimcount>1)
        error('getMot:extraXLDdimEl2','extraXLDdim may not contain the same dimension number twice');
    elseif size(extraXLDdim,2) == 1
        extraXLDdim = extraXLDdim';
    end
    
    if any(diff(extraXLDdim)<0) % User didn't specify dimensions in order
        [extraXLDdim,perm] = sort(extraXLDdim,'ascend');
        extraXLD = extraXLD(perm);
    end
    
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

% Check if sync indices exist for MarkerData and ForceData
if ~isfield(C3Ddata.Marker,'MarkerSyncIdx') || ~isfield(C3Ddata.Force,'ForceSyncIdx')
    warning('getMot:syncIdx',['No sync indices found in' infilename '. Assuming synchronized marker-force data']);
    
    markerSyncIdx = [1 size(C3Ddata.Marker.MarkerData,1)];
    forceSyncIdx = [1 size(C3Ddata.Force.ForceData,1)];
else
    markerSyncIdx = C3Ddata.Marker.MarkerSyncIdx;
    forceSyncIdx = C3Ddata.Force.ForceSyncIdx;
end

% Check if sync indices exist for OtherData
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
% columnNames = {'time' ...
%     'l_ground_force_vx' 'l_ground_force_vy' 'l_ground_force_vz' ... % Left GRF Vector in specific body CRF
%     'l_ground_torque_x' 'l_ground_torque_y' 'l_ground_torque_z' ... % Left Moments
%     'ground_force_vx' 'ground_force_vy' 'ground_force_vz' ... % Right GRF Vector in specific body CRF
%     'ground_torque_x' 'ground_torque_y' 'ground_torque_z' ... % Right Moments
%     'l_ground_force_px' 'l_ground_force_py' 'l_ground_force_pz' ... % Left COP
%     'ground_force_px' 'ground_force_py' 'ground_force_pz' }; % Right COP

%like in template .mot for GRF (Force COP l_Force l_COP Moment l_Moment)
columnNames = {'time' ...
    'l_ground_force_vx' 'l_ground_force_vy' 'l_ground_force_vz' ... % Left GRF Vector in specific body CRF
    'l_ground_force_px' 'l_ground_force_py' 'l_ground_force_pz' ... % Left COP
    'ground_force_vx' 'ground_force_vy' 'ground_force_vz' ... % Right GRF Vector in specific body CRF
    'ground_force_px' 'ground_force_py' 'ground_force_pz' ... % Right COP
    'l_ground_torque_x' 'l_ground_torque_y' 'l_ground_torque_z' ... % Left Moments
     'ground_torque_x' 'ground_torque_y' 'ground_torque_z'}; % Right Moments

% Add additional external loads, if any
if ~isempty(extraXLD)
    
    if isempty(extraXLDdim)
        for ixld = 1:length(extraXLD);
            columnNames(end+1:end+3) = {[extraXLD{ixld} '_x'] , [extraXLD{ixld} '_y'] , [extraXLD{ixld} '_z']};
        end
    else
        
        n = 1; 
        nInGroup = zeros(1,numel(extraXLDdim));
        while any((extraXLDdim-3*(n-1))>0)
            nInGroup(n) = sum((extraXLDdim - 3*n)<=0) - sum((extraXLDdim - 3*(n-1))<=0);
            n = n + 1;
        end
        nInGroup(nInGroup<=0) = [];
        
        for ixld = cumsum([0 nInGroup(1:end-1)])+1 % First names of groups
            columnNames(end+1:end+3) = {[extraXLD{ixld} '_x'] , [extraXLD{ixld} '_y'] , [extraXLD{ixld} '_z']};
        end
        
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

    if isempty(extraXLDdim)
        otherData = zeros(length(otherIdx),3*length(addIdx));
        otherData(:,1 + (0:3:3*length(addIdx)-1) ) = C3Ddata.Other.OtherData(otherIdx,addIdx); % Each element in extraXLD gets its own xyz input, with yz being 0
    else
        otherData = zeros(length(otherIdx),3*length(nInGroup));
        otherData(:,extraXLDdim) = C3Ddata.Other.OtherData(otherIdx,addIdx);
    end
    
    % Apply permVec
    osiz = size(otherData,2);
    permadd = repmat(zeros(1,3),[1 osiz/3]);
    permadd(1:3:osiz) = 0:3:3*(osiz/3 - 1);
    permadd = cumsum(permadd);
    permall = repmat(permvec,[1 osiz/3]) + permadd;
    otherData = otherData(:,permall);
    
else
    otherData = [];
end

% COP data (zeros, as we apply a force AND a torque, 
% unless the forceplate CRF has an offset wrt the MoCap CRF)
% TODO: allow for time varying plate offsets (timeseries in structure)
%C3Ddata.Force.ForcePlateOffset=1; %Include COP
%if ~isfield(C3Ddata.Force,'ForcePlateOffset')
%    copData = zeros(size(forceData,1),6);
%else
copData = zeros(size(forceData,1),6);
copStruct = getCOP(forceData,fthresh);
copData(:,1:3) = copData(:,1:3)+copStruct.COPL;
copData(:,4:6) = copData(:,4:6)+copStruct.COPR;
%end

%Zeros for all M but My (free torque, it has the vertical (Z) (Y in OpenSim) component only: http://www.kwon3d.com/theory/grf/grf.html)
forceData(:,4) = zeros(size(forceData,1),1); %MxL
forceData(:,6) = zeros(size(forceData,1),1); %MzL
forceData(:,10) = zeros(size(forceData,1),1); %MxR
forceData(:,12) = zeros(size(forceData,1),1); %MzR



% Add time column
writeData = [(0:size(forceData,1)-1)'./markerFrameRate forceData copData otherData];


%like in template .mot for GRF (Force COP l_Force l_COP Moment l_Moment)
writeData = [writeData(:,1:4) writeData(:,14:16) writeData(:,8:10) writeData(:,17:19) writeData(:,5:7) writeData(:,11:13)];

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
    'inDegrees=yes\n' ...
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