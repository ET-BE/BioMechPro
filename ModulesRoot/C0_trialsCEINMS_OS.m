function [Datastr] = C0_trialsCEINMS_OS(Datastr,ceinmsFolder,osFolder,emgFolder)
% gBMPDynUI ceinmsFolder=1; osFolder=1; emgFolder=1;
% 
% For each trial, the function creates an .xml file containing the paths that CEINMS needs when running calibration and execution (EMG,
% ID, MA, etc)
% 
% This is for trials using the Optical data (joint angles were obtained via OpenSim (OS))
%% Check

% Get data from .Info field (assumed there)
rootfolder = Datastr.Info.SubjRoot;
bckslsh = strfind(rootfolder,'\');
if isempty(bckslsh)
    bckslsh = 0;
end
trial = Datastr.Info.Trial;
savename = [rootfolder(bckslsh(end)+1:end) trial];
osmastofolder = [rootfolder '\' osFolder '\' savename 'MA'];
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
osmaFilePath = [rootfolder '\' osFolder '\' savename 'MA'];
osidstoFilePath = [rootfolder '\' osFolder '\' savename 'ID.sto']; % Path to ID output (.sto file)
trialsGenSetPath = [rootfolder '\' ceinmsFolder '\trials\' trialsGenSet]; % 
emgFilePath = [rootfolder '\' emgFolder '\' savename '\' 'emg.mot']; % IK output


%% CEINMS

% create the xlm files of trials
getCEINMStrials(osmaFilePath,osidstoFilePath,trialsGenSetPath,emgFilePath);



end