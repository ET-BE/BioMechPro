function [Datastr] = OS6_OSMA_IMU(Datastr,osInstallPath,osFolder,imuFolder,permvec)
% gBMPDynUI osInstallPath=1; osFolder=1; imuFolder=1; permVec=1;
%
% INPUT)
% - Datastr: with no specific fields required
% 
% - osInstallPath: string, full path of the OpenSim directory. Example:
% 'C:\Program Files\OpenSim\OpenSim 3.3';
% 
% - osFolder, string, specifying the folder in which the OpenSim files will
% be stored, relative to the subject root folder. Example:
% 'OS'
% 
% - permVec, vector, containing 3 elements to permute the data dimensions 
% (xyz) from OpenSim to whatever order you're using in the structure. 
% When left empty, no permutation occurs.
% Default OpenSim: x = walking direction, z = to the right, y = upward
% Example:
% [2 3 1];
% 
% OUTPUT)
% - Datastr: structure, with added fields:
% - MA files

% NOTES)
% Same as OS6_OSMA but for the inertial data
%
%% Check
if ~isfield(Datastr,'IMU') % if no IMU data
    warning(['No IMU data in trial ' Datastr.Info.Trial '. Skipping.']);
    return
end


%% Get info

% Get data from .Info field (assumed there)
rootfolder = Datastr.Info.SubjRoot;
% bckslsh = strfind(rootfolder,'\');
% if isempty(bckslsh)
%     bckslsh = 0;
% end

trial = Datastr.Info.Trial;


% Get all .mot filenames in IMU folder
imufiles = dir([rootfolder '\' imuFolder '\*.mot']);
imufilenames = cell(1,length(imufiles));
imufilenamescmp = cell(1,length(imufiles));
for ifile = 1:length(imufiles)
    imufilenames{ifile} = imufiles(ifile).name;
    
    % Reverse string and remove '.mat', for comparison with strncmpi
    imufilenamescmp{ifile} = imufiles(ifile).name(end-4:-1:1);
end

% Check if storage filenames specified by Trials field exist.
% If not, skip the thing
match = strncmpi( trial(end:-1:1),imufilenamescmp,length(trial) );

if sum(match)==1
    subjimufile = imufilenames{match};
else
    warning('OS6_OSMA_IMU:nofile',['Could not find unique imu file ' trial '. Skipping.']);
    return;
end

% Get data from .Info field (assumed there)
rootfolder = Datastr.Info.SubjRoot;
bckslsh = strfind(rootfolder,'\');
if isempty(bckslsh)
    bckslsh = 0;
end

trial = Datastr.Info.Trial;
savename = [rootfolder(bckslsh(end)+1:end) trial];

osMod = Datastr.Info.subjosmodfile;
maGenSet = Datastr.Info.subjosmasetfile;

% Create paths
osModPath = [rootfolder '\' osFolder '\' osMod]; % Subject model
maGenSetPath = [rootfolder '\' osFolder '\' maGenSet]; % MA general settings

ikFilePath = [rootfolder '\' imuFolder '\' subjimufile]; % IK input


%% Do MA

% Do MA
getOSMA_IMU(osInstallPath,osModPath,maGenSetPath,ikFilePath);

% Get path to created .sto file
osmastofolder = [rootfolder '\' imuFolder '\' savename 'MA']; % Path to MA output (.sto file)
% Store in Datastr structure
Datastr = getOSMAinM(Datastr,osmastofolder,permvec);


end