function [Datastr] = S10_toOS(Datastr, osFolder, permVec, removeEndFrames)
% gBMPDynUI osFolder=1; permVec=1; removeEndFrames=1;
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
% beginEndRemove: remove 90 frames at the beginning and at the end, to
% remove the artifacts of sychronization
%
% OUTPUT)
% - None
% 
% NOTES)

%% Get info

subjroot = Datastr.Info.SubjRoot;
trial = Datastr.Info.Trial;

if isempty(permVec)
    permVec = 1:3;
end

%% remove fields if exist
if isfield(Datastr.Resample, 'Sych')
    if isfield(Datastr.Resample.Sych, 'IKAngData')
        Datastr.Resample.Sych = rmfield(Datastr.Resample.Sych, 'IKAngData');
    end
    
    if isfield(Datastr.Resample.Sych, 'ForcePlateGRFData')
        Datastr.Resample.Sych = rmfield(Datastr.Resample.Sych, 'ForcePlateGRFData');
    end
    
end

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

% Create .sto file (kinematic data)
if isfield(Datastr.Resample, 'Marker')
    getTrc(Datastr,[subjroot '\' osFolder '\' subjroot(bckslsh(end)+1:end) trial], permVec, removeEndFrames)
else
    warning(['No field Marker in trial' trial '. Skipping writing to .trc file. Unable to do IK.']);
end

if isfield(Datastr.Resample,'Force')  % Create .mot file (kinetic data)
    Datastr = getMot(Datastr,[subjroot '\' osFolder '\' subjroot(bckslsh(end)+1:end) trial], removeEndFrames);
else
    warning(['No field Force in trial' trial '. Skipping writing to .mot file. Unable to do ID.']);
end

% Set empty return value (to prevent saving by UI)
% Datastr = [];
end