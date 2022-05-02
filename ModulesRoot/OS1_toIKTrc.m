function [Datastr] = OS1_toIKTrc(Datastr,osFolder,permVec,expFlag)
% gBMPDynUI osFolder=1; permVec=1; expFlag=1;
% 
% Create .mot files and .trc files for use in OpenSim
% 
% INPUT)
% - Datastr, the data structure with at least the fields:
% (TODO)
% 
% - osFolder: string, specifying the folder in which the OpenSim files will
% be stored, relative to the subject root folder. Example:
% 'OS'
% 
% - permVec: vector, containing 3 elements to permute the data dimensions 
% (xyz) to OpenSim coordinates. When left empty, no permutation occurs.
% Default OpenSim: x = walking direction, z = to the right, y = upward
% Example:
% [2 3 1];
% 
% - expFlag: scalar, flag to specify which data is exported.
% 1) Marker data
% 2) Probe data (which will then be treated as marker data)
% 3) Both
% Default is 1
% 
% 
% OUTPUT)
% - None
% 
% NOTES)

%% Get info

subjroot = Datastr.Info.SubjRoot;
trial = Datastr.Info.Trial;

if ~isfield(Datastr,'Marker') % if no VZ data
    warning(['No VZ data in file ' Datastr.Info.Trial '. Skipping.']);
    return
end

if isempty(permVec)
    permVec = 1:3;
end
if isempty(expFlag)
    expFlag = 1;
end

%% Folder 

% Check if folder exist, if not create new
if exist([subjroot '\' osFolder],'dir') ~= 7
    mkdir(subjroot,osFolder);
end

%% Export files

bckslsh = strfind(subjroot,'\');
if isempty(bckslsh)
    bckslsh = 0;
end

% Create .trc file (kinematic data)
getTrc(Datastr,[subjroot '\' osFolder '\' subjroot(bckslsh(end)+1:end) trial],expFlag,permVec)

% Set empty return value (to prevent saving by UI)
Datastr = [];

end