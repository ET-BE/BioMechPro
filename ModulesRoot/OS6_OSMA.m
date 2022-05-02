function [Datastr] = OS6_OSMA(Datastr,osInstallPath,osFolder,permvec)
% gBMPDynUI osInstallPath=1; osFolder=1; permVec=1;
% 
% Run muscle analysis (MA) using OpenSim, and port the data to the 
% structure
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
% 
%% Check
%Rough method to check that we should have .mot files with joint angles
%after IK. It doesn't check if we performed IK!
if ~isfield(Datastr,'Marker') % if no VZ data
    warning(['No .trc in trial ' Datastr.Info.Trial '. Skipping.']);
    return
end

if ~isfield(Datastr,'Force') % if no FP data
    warning(['No XLD.mot (loads) in trial ' Datastr.Info.Trial '. Skipping.']);
    return
end

%% Get info

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

ikFilePath = [rootfolder '\' osFolder '\' savename 'IK.mot']; % IK input


%% Do MA

% Do MA 
getOSMA(osInstallPath,osModPath,maGenSetPath,ikFilePath);

% Get path to created MA folder
osmastofolder = [rootfolder '\' osFolder '\' savename 'MA']; % Path to MA output (.sto file)
% Store in structure
Datastr = getOSMAinM(Datastr,osmastofolder,permvec);


end