function [Datastr] = S13_toOS_Portable(Datastr, osFolder, osimDoFs)
% gBMPDynUI osFolder=1; osimDoFs=1;
% 
% Create .sto files and .mot files for use in OpenSim, from the portable
% system.
% 
% INPUT)
% - Datastr, the data structure with at least the fields:
% (TODO)
% 
% - osFolder: string, specifying the folder in which the OpenSim files will
% be stored, relative to the subject root folder. Example:
% 'OS'
% 
% - Sign: scalar, whether the joint defintion is the same between OpenSim
% and IMU systems
%
% OUTPUT)
% - None
% 
% NOTES)
% - PortableFlag: string, flag to specify if the portable data will be
% exported. Default is 'false'

%% check input
if strcmp(num2str(osimDoFs), '23')
    IKheaders = ["time",	"pelvis_tilt",	"pelvis_list",	"pelvis_rotation",...
               "pelvis_tx",	"pelvis_ty",	"pelvis_tz",	"hip_flexion_r",...
               "hip_adduction_r", "hip_rotation_r", "knee_angle_r", "ankle_angle_r",...
               "subtalar_angle_r", "mtp_angle_r", "hip_flexion_l", "hip_adduction_l",...
               "hip_rotation_l", "knee_angle_l", "ankle_angle_l", "subtalar_angle_l",...
               "mtp_angle_l",	"lumbar_extension",	"lumbar_bending",	"lumbar_rotation"];
else
    error('new headers are needed, if the Osim model is not Gait2392')
end

%% remove fields if exist
if isfield(Datastr.Resample, 'Sych')
    if isfield(Datastr.Resample.Sych, 'IMUAngData')
        Datastr.Resample.Sych = rmfield(Datastr.Resample.Sych, 'IMUAngData');
    end
    
    if isfield(Datastr.Resample.Sych, 'InsoleGRFData')
        Datastr.Resample.Sych = rmfield(Datastr.Resample.Sych, 'InsoleGRFData');
    end
end

%% extract subject information

%% Get info

subjroot = Datastr.Info.SubjRoot;
trial = Datastr.Info.Trial;

%% Folder 

% Check if folder exist, if not create new
if ~ exist([subjroot '\' osFolder],'dir')
    mkdir(subjroot,osFolder);
end

%% Export files
bckslsh = strfind(subjroot,'\');
if isempty(bckslsh)
    bckslsh = 0;
end

if isfield(Datastr.Resample, 'IMU')    
    Datastr = getIMUMot(Datastr,[subjroot '\' osFolder '\DataFiles\' subjroot(bckslsh(end)+1:end) trial], IKheaders);
else
    warning(['No field IMU in trial' trial '. Skipping writing to .sto file. Unable to do portable ID.']);
end

% Create .mot file (kinetic data)
if isfield(Datastr.Resample,'Insole')
    Datastr = getInsoleMot(Datastr,[subjroot '\' osFolder '\DataFiles\' subjroot(bckslsh(end)+1:end) trial]);
else
    warning(['No field Insole in trial' trial '. Skipping writing to .mot file. Unable to do portable ID.']);
end

% Set empty return value (to prevent saving by UI)
% Datastr = [];
end