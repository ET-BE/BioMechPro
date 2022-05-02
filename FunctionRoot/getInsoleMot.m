function [C3Ddata] = getInsoleMot(C3Ddata, filename)
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
% sychTime: the delta T of the force schronization. Insole data is deltaT
% time behand the treadmill data
% 
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

% Data column header
columnNames = {'time' ...
    'l_ground_force_vx' 'l_ground_force_vy' 'l_ground_force_vz' ... % Left GRF Vector in specific body CRF
    'l_ground_force_px' 'l_ground_force_py' 'l_ground_force_pz' ... % Left COP
    'ground_force_vx' 'ground_force_vy' 'ground_force_vz' ... % Right GRF Vector in specific body CRF
     'ground_force_px' 'ground_force_py' 'ground_force_pz' ... % Right COP   
    'l_ground_torque_x' 'l_ground_torque_y' 'l_ground_torque_z' ... % Left Moments
    'ground_torque_x' 'ground_torque_y' 'ground_torque_z' }; % Right Moments

nRowsIKAng = size(C3Ddata.Resample.Sych.IKAngData, 1);  % make sure the insole forces have the same rows as the IKAng

if isfield(C3Ddata.Resample.Sych, 'DeltaT')
    sychT = C3Ddata.Resample.Sych.DeltaT.Force;
else
    sychT = 0.1;
end

left_force = interp1((0:size(C3Ddata.Resample.Insole.Left, 1)-1)/FrameRate,...
             C3Ddata.Resample.Insole.Left(:, 23:25),...
             (0:nRowsIKAng-1)/FrameRate + sychT,...
             'linear', 'extrap');
         
right_force = interp1((0:size(C3Ddata.Resample.Insole.Right, 1)-1)/FrameRate,...
             C3Ddata.Resample.Insole.Right(:, 23:25),...
             (0:nRowsIKAng-1)/FrameRate + sychT,...
             'linear', 'extrap');

%% Get number of force samples equal to number of marker samples


%% Collect FORCE, MOMENT and COP data to write

% Collect force and moment data
% permall = [permvec permvec+3 permvec+6 permvec+9];  % Assumed 12 channel MGRF
forceData = zeros(nRowsIKAng, 12);
forceData(:, 2) = left_force(:, 1);
forceData(:, 8) = right_force(:, 1);

% COP data 
% unless the forceplate CRF has an offset wrt the MoCap CRF)
% TODO: allow for time varying plate offsets (timeseries in structure)
if C3Ddata.Info.subjShoesize == 42 || C3Ddata.Info.subjShoesize == 43
    insole_size  = [274.2, 97.5]/1000;
elseif C3Ddata.Info.subjShoesize == 38 || C3Ddata.Info.subjShoesize == 39
    insole_size  = [248.6, 90.2]/1000;
end

copData = zeros(nRowsIKAng, 6);
% relocate the origin of the cop in the x direction to the rear of the
% insole, and change the right y cop direction to the oppsite, to the right
% side of the body.
copData(:, [1, 3]) = [(left_force(:, 2)+0.5).*insole_size(1),...
                       left_force(:, 3).*insole_size(2)];
copData(:, [4, 6]) = [(right_force(:, 2)+0.5).*insole_size(1),...
                      -right_force(:, 3).*insole_size(2)];  

% Add time column
writeData = [(0:size(forceData,1)-1)'./FrameRate forceData(:, 1:3) copData(:,1:3) forceData(:, 7:9) copData(:,4:6) forceData(:, 4:6) forceData(:, 10:12)];

%% Create file and write header
fid = fopen([filename 'XLD_INSOLE.mot'],'w');
if fid == -1
    error('getMot:FileID',['Cannot open ' infilename '.mot for writing. It might be in use.']);
end

% General header (SIMM header only)
fprintf(fid, [infilename '\n' ...
    'version=1\n' ...
	'nRows=' num2str(nRowsIKAng) '\n' ...
    'nColumns=' num2str(length(columnNames)) '\n' ...
    'inDegrees=no\n' ...
    'endheader\n']);

% Column headers
for iCol = 1:length(columnNames)
    fprintf(fid, [columnNames{iCol} '\t']);
end
fprintf(fid, '\n');

%% Write data
writeStr = regexprep(mat2str(writeData),{'[',']',' ',';'},{'','','\t','\n'});
fprintf(fid,writeStr);

%% Clean up

fclose(fid);

disp([infilename 'XLD_INSOLE.mot created in ' pathname]);

C3Ddata.Resample.Sych.InsoleGRFData = writeData;
C3Ddata.Resample.Sych.InsoleGRFDataLabel = columnNames;

end