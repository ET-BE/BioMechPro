function [Datastr] = C0_trialsCEINMS_IMU(Datastr,ceinmsFolder,imuFolder,emgFolder)
% gBMPDynUI ceinmsFolder=1; imuFolder=1; emgFolder=1;
% 
% For each trial, the function creates an .xml file containing the paths that CEINMS needs when running calibration and execution (EMG,
% ID, MA, etc)
% 
% This is for trials using the inertial data (joint angles were obtained via IMU)
%% Check

% Get data from .Info field (assumed there)
rootfolder = Datastr.Info.SubjRoot;
bckslsh = strfind(rootfolder,'\');
if isempty(bckslsh)
    bckslsh = 0;
end
trial = Datastr.Info.Trial;
savename = [rootfolder(bckslsh(end)+1:end) trial];
osmastofolder = [rootfolder '\' imuFolder '\' savename 'MA'];
emgtofolder = [rootfolder '\' emgFolder '\' savename];

if 7~=exist([osmastofolder ],'dir') % if no MA data
    warning(['No Muscle Analysis performed in trial ' Datastr.Info.Trial '. Skipping.']);
    return
end
if 7~=exist([emgtofolder],'dir') % if no EMG data
    warning(['No EMG normalization in trial ' Datastr.Info.Trial '. Skipping.']);
    return
end

trialsGenSet = Datastr.Info.subjceinmstrialssetfile;

% Create paths
osmaFilePath = [rootfolder '\' imuFolder '\' savename 'MA'];
idstoFilePath = [rootfolder '\' ceinmsFolder '\trials\' 'templateID.sto']; % Path to ID output (.sto file), non existent with IMU but needed for CEINMS not to crash
trialsGenSetPath = [rootfolder '\' ceinmsFolder '\trials\' trialsGenSet]; % 
emgFilePath = [rootfolder '\' emgFolder '\' savename '\' 'emg.mot']; % IK output


%% CEINMS

% create the .xlm files of trials
getCEINMStrials(osmaFilePath,idstoFilePath,trialsGenSetPath,emgFilePath);

end