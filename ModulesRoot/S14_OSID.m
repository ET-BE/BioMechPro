function [Datastr] = S14_OSID(Datastr,osInstallPath,osFolder,permvec)
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
if ~isfield(Datastr,'Marker') % if no VZ data
    warning(['No .trc in trial ' Datastr.Info.Trial '. Skipping.']);
    return
end

if ~isfield(Datastr,'Force') % if no VZ data
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
idGenSet = Datastr.Info.subjosidsetfile;
xldGenSet = Datastr.Info.subjosidxldfile;

% Create paths
osModPath = [rootfolder '\' osFolder '\' osMod]; % Subject model
idGenSetPath = [rootfolder '\' osFolder '\' idGenSet]; % ID general settings
xldGenSetPath = [rootfolder '\' osFolder '\' xldGenSet]; % External load ID settings

ikFilePath = [rootfolder '\' osFolder '\DataFiles\' savename 'IK.mot']; % IK output
xldFilePath = [rootfolder '\' osFolder '\DataFiles\' savename 'XLD.mot']; % External load file

%% Do ID

% Do ID
getOSID(osInstallPath,osModPath,idGenSetPath,xldGenSetPath,ikFilePath,xldFilePath);

% Get path to created ID.sto file
osidsto = [rootfolder '\' osFolder '\DataFiles\' savename 'ID.sto']; % Path to ID output (.sto file)

% Store in C3Ddata structure
Datastr = getOSIDinM(Datastr,osidsto,permvec);



end