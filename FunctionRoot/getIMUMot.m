function [C3Ddata] = getIMUMot(C3Ddata, filename, IKheaders)
%% getSto
% Get .mot joint angle file from IMU data, for use in portable ID in OpenSim
%
% INPUT)
% C3Ddata : Data structure containing at least the fields:
% [...]
% 
% trialname : Path and name of the .sto file to be created.
% If no path is provided the file is saved in the current folder.
%
% sychTime: the delta T of the force schronization. Insole data is deltaT
% time behand the treadmill data
% 
% OUTPUT)
% No direct output
% A .sto file is created in the destination provided by trialname

% Mark Vlutters - July 2015 - Enschede
% Modify by Huawei Wang - July 2021 - Enschede

%% Check input

% Check trialname
if ~ischar(filename)
    error('getSto:trialname','Input trialname must be a string');
end

%% Collect header to write
if isempty(IKheaders)
    IKheaders = ["time",	"pelvis_tilt",	"pelvis_list",	"pelvis_rotation",...
               "pelvis_tx",	"pelvis_ty",	"pelvis_tz",	"hip_flexion_r",...
               "hip_adduction_r", "hip_rotation_r", "knee_angle_r", "ankle_angle_r",...
               "subtalar_angle_r", "mtp_angle_r", "hip_flexion_l", "hip_adduction_l",...
               "hip_rotation_l", "knee_angle_l", "ankle_angle_l", "subtalar_angle_l",...
               "mtp_angle_l",	"lumbar_extension",	"lumbar_bending",	"lumbar_rotation"];
end
           
IMUheaders = C3Ddata.IMU.IMUDataLabel;
% {'pelvis_tilt','pelvis_list','pelvis_rotation','pelvis_tx','pelvis_ty','pelvis_tz',...
% 'hip_flexion_r','hip_adduction_r','hip_rotation_r','knee_angle_r','knee_adduction_r',...
% 'knee_rotation_r','ankle_angle_r','ankle_adduction_r','subtalar_angle_r','mtp_angle_r',...
% 'hip_flexion_l','hip_adduction_l','hip_rotation_l','knee_angle_l','knee_adduction_l',...
% 'knee_rotation_l','ankle_angle_l','ankle_adduction_l','subtalar_angle_l','mtp_angle_l',...
% 'lumbar_extension','lumbar_bending','lumbar_rotation'}
index_select = [];
for dof = IKheaders(2:end)
    index_select = [index_select, find(strcmp(IMUheaders, dof))];
end

writeHeader = IKheaders;

%% Collect some info
FrameRate = C3Ddata.Resample.FrameRate;
% IMUFrameRate = C3Ddata.IMU.IMUFrameRate;
% nDoFs = length(writeHeader) - 1;

if isfield(C3Ddata.Resample.Sych, 'DeltaT')
    sychT = C3Ddata.Resample.Sych.DeltaT.Motion;
else
    sychT = 0.05;
end


nFramesIKAng = size(C3Ddata.Resample.Sych.IKAngData, 1);
nFramesIMU = size(C3Ddata.Resample.IMU, 1);
% IMUnFrames = size(C3Ddata.IMU.IMUData, 1);
% sychronize with marker data

IMUData = interp1((0:nFramesIMU-1)/FrameRate, C3Ddata.Resample.IMU(1:nFramesIMU, index_select),...
          (0:nFramesIKAng-1)/FrameRate + sychT, 'linear', 'extrap');
CoMData = interp1((0:nFramesIMU-1)/FrameRate, C3Ddata.Resample.CoM(1:nFramesIMU, :),...
          (0:nFramesIKAng-1)/FrameRate + sychT, 'linear', 'extrap');

% If full path is supplied, take last part for name inside trc file
if ~isempty(strfind(filename, '\'))
    foo = strfind(filename, '\');
    infilename = filename(foo(end)+1:end);
    pathname = filename(1:foo(end));
else
    infilename = filename;
    pathname = pwd;
end

%% Collect data to write

msiz = size(IMUData);
writeData = [...
            (0:msiz(1)-1)'./FrameRate , ...
             IMUData];

%% Create file and headers
% File will contain all markers and probe positions
% Not all might have to be associated with the scaling of the model

% File
fid = fopen([filename 'IK_IMU.mot'],'w'); % note: w also discards existing content, if any
if fid == -1
    error('getSto:FileID',['Cannot open ' infilename '.sto for writing. It might be in use.']);
end

% General header
fprintf(fid, [infilename 'Coordinates\n' ...
    'version=1\n' ...
	'nRows=' num2str(msiz(1)) '\n' ...
    'nColumns=' num2str(msiz(2)+1) '\n' ...
    'inDegrees=yes\n' ...
    '\n' ...
    'Units are S.I. units (second, meters, Newtons, ...)\n' ...
    'Angles are in degrees.\n\n' ...
    'endheader\n']);

% Column headers
for iCol = 1:msiz(2)+1
    fprintf(fid, [IKheaders{iCol} '\t']);
end
fprintf(fid, '\n');

%% Write data

writeStr = regexprep(mat2str(writeData),{'[',']',' ',';'},{'','','\t','\n'});
fprintf(fid,writeStr);

%% Clean up

fclose(fid);

disp([infilename 'IK_IMU.sto created in ' pathname]);

C3Ddata.Resample.Sych.IMUAngData = writeData;
C3Ddata.Resample.Sych.IMUAngDataLabel = IKheaders;
C3Ddata.Resample.Sych.CoM = CoMData;

end