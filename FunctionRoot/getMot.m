function [C3Ddata] = getMot(C3Ddata, filename, removeEndFrames)
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

% Frame rates
FrameRate = C3Ddata.Resample.FrameRate;
% forceFrameRate = C3Ddata.Force.ForceFrameRate;

% Data column header

columnNames = {'time' ...
    'l_ground_force_vx' 'l_ground_force_vy' 'l_ground_force_vz' ... % Left GRF Vector in specific body CRF
    'l_ground_force_px' 'l_ground_force_py' 'l_ground_force_pz' ... % Left COP
    'ground_force_vx' 'ground_force_vy' 'ground_force_vz' ... % Right GRF Vector in specific body CRF
     'ground_force_px' 'ground_force_py' 'ground_force_pz' ... % Right COP   
    'l_ground_torque_x' 'l_ground_torque_y' 'l_ground_torque_z' ... % Left Moments
    'ground_torque_x' 'ground_torque_y' 'ground_torque_z' }; % Right Moments

% Data size
nRows = size(C3Ddata.Resample.Marker, 1) - removeEndFrames;  % make sure the treadmill forces has the same frames as markers
nColumns = length(columnNames);
%% Collect FORCE, MOMENT and COP data to write

% Collect force and moment data
% permall = [permvec permvec+3 permvec+6 permvec+9];  % Assumed 12 channel MGRF
% forceData = C3Ddata.Force.ForceData(forceIdx,permall);
% momentLData = C3Ddata.Resample.Force.Left(:,[4:6]);
% forceLData = C3Ddata.Resample.Force.Left(:,[1:3]);
% momentRData = C3Ddata.Resample.Force.Right(:,[4:6]);
% forceRData = C3Ddata.Resample.Force.Right(:,[1:3]);
forceData = [C3Ddata.Resample.Force.Left(1:nRows, 1:6) C3Ddata.Resample.Force.Right(1:nRows, 1:6)];

% COP data 
% unless the forceplate CRF has an offset wrt the MoCap CRF)
% TODO: allow for time varying plate offsets (timeseries in structure)
if ~isfield(C3Ddata.Force,'ForcePlateOffset')
    copData = zeros(nRows,6);

else
    copData = repmat(C3Ddata.Force.ForcePlateOffset,[nRows 1]);
end
copStruct = getCOP(forceData,fthresh);  % the recorded CoP is not at the same coordiante with Force/Moment, especially in X direction
copData(:,1:3) = copData(:,1:3)+copStruct.COPL;
copData(:,4:6) = copData(:,4:6)+copStruct.COPR;

%Zeros for all M but My (free torque, it has the vertical (Z) (Y in OpenSim) component only: http://www.kwon3d.com/theory/grf/grf.html)
forceData(:,4) = zeros(size(forceData,1),1); %MxL
forceData(:,6) = zeros(size(forceData,1),1); %MzL
forceData(:,10) = zeros(size(forceData,1),1); %MxR
forceData(:,12) = zeros(size(forceData,1),1); %MzR

% if the CoP is NaN, make all the forces and moments equal zero
forceData(isnan(copData(:, 1)), 1:6) = 0;
forceData(isnan(copData(:,4)), 7:12) = 0;

% resign CoP to zeros from NaNs.
copData(isnan(copData)) =0;

% Add time column
writeData = [(0:nRows-1)'./FrameRate forceData(:, 1:3) copData(:,1:3) forceData(:, 7:9) copData(:,4:6) forceData(:, 4:6) forceData(:, 10:12)];

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

C3Ddata.Resample.Sych.ForcePlateGRFData = writeData;
C3Ddata.Resample.Sych.ForcePlateGRFDataLabel = columnNames;

end