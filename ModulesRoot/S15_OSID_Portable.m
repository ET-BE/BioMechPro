function [Datastr] = S15_OSID_Portable(Datastr,osInstallPath,osFolder,permvec)
% gBMPDynUI osInstallPath=1; osFolder=1; permVec=1;
% 
% INPUT)
% 
% OUTPUT)
% 
% NOTES)
% 
% 
%% Check
if ~isfield(Datastr,'IMU') % if no VZ data
    warning(['No _IMU.mot in trial ' Datastr.Info.Trial '. Skipping.']);
    return
end

if ~isfield(Datastr,'Insole') % if no VZ data
    warning(['No XLD_INSOLE.mot (loads) in trial ' Datastr.Info.Trial '. Skipping.']);
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
idGenSet = Datastr.Info.subjosidsetfile;
xldGenSet = Datastr.Info.subjosidxldfile_p;

% Create paths
osModPath = [rootfolder '\' osFolder '\' osMod]; % Subject model
idGenSetPath = [rootfolder '\' osFolder '\' idGenSet]; % ID general settings
xldGenSetPath = [rootfolder '\' osFolder '\' xldGenSet]; % External load ID settings

ikFilePath = [rootfolder '\' osFolder '\DataFiles\' savename 'IK_IMU.mot']; % IK output
xldFilePath = [rootfolder '\' osFolder '\DataFiles\' savename 'XLD_INSOLE.mot']; % External load file

%% Do ID

% Do ID
getOSID_Portable(osInstallPath,osModPath,idGenSetPath,xldGenSetPath,ikFilePath,xldFilePath);

% Get path to created ID.sto file
osidsto = [rootfolder '\' osFolder '\DataFiles\' savename 'ID_Portable.sto']; % Path to ID output (.sto file)

% Store in C3Ddata structure
Datastr = getOSIDinM_Portable(Datastr,osidsto,permvec);

end