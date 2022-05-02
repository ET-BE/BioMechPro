function [Datastr] = OS1_toOS(Datastr,osFolder,permVec,expFlag,extraXLDforce,extraXLDposition,extraXLDdim,extraNames)
% gBMPDynUI osFolder=1; permVec=1; expFlag=1; extraXLD=1; extraXLDposition =1; extraXLDdim=1; extraNames=1;
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
% - extraXLD: string or cell of strings, containing label names in
% OtherDataLabel. Corresponding columns of OtherData will also be written
% to the file.
% 
% - extraXLDdim : numerical vector, containing a dimension index 1,2 or 3
% specifying in which OpenSim dimension the extraXLD channels are stored.
% The numbers 1-3 give suffix 'x','y', and 'z' respectively.
% The two other not specified dimensions will contain only zeros.
% The number of elements in extraXLDg must be the same as the number of 
% channels specified by extraXLD. If extraXLDdim is not specified, it will
% default to 1 for every channel in extraXLD.
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

if isfield(Datastr,'Force')
    
    % Do some data manipulation for the extra load
%     if ~isempty(extraXLD) % Convert Nm to N and flip sign
%         for iExtra = 1:length(extraXLD)
%             idx(iExtra) = find(strcmp(Datastr.Other.OtherDataLabel,extraXLD(iExtra)));
%         end
    
    % Create .mot file (kinetic data)
        getMot(Datastr,[subjroot '\' osFolder '\' subjroot(bckslsh(end)+1:end) trial],permVec);
    
else
    warning(['No field Force in trial' trial '. Skipping writing to .mot file. Unable to do ID.']);
    return;
end

% Set empty return value (to prevent saving by UI)
Datastr = [];

end