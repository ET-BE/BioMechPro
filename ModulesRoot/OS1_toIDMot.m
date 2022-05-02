function [Datastr] = OS1_toIDMot(Datastr,osFolder,permVec,extraXLD,extraXLDdim)
% gBMPDynUI osFolder=1; permVec=1; extraXLD=1; extraXLDdim=1;
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
% permvec : permutation vector to put data in OpenSim coordinates.
% (x: walking direction, z: to the right, y: upward)
% Default is [1 2 3] (no permutation)
% This vector is applied to both the ForceData in the data structure, as
% well as the extraXLD that is optionally specified.
% 
% extraXLD : string or cell of strings, containing label names in
% OtherDataLabel. Corresponding columns of OtherData will also be written
% to the file as external loads.
% 
% extraXLDdim : numerical vector, containing one or more scalar dimension 
% indices, specifying in which OpenSim dimension the extraXLD channels are 
% stored. The numbers 1:3 give suffix 'x','y', and 'z' respectively.
% Use 4:6 for the next group of variables, 7:9 for the next, etc.
% Not specified dimensions will contains zeroes.
% The number of elements in extraXLDdim must be the same as the number of 
% channels specified by extraXLD. If extraXLDdim is not specified, each
% channel specified by extraXLD will be stored in the x dimension.
% 
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

%% Check

if ~isfield(Datastr,'Force') % if no VZ data
    warning(['No FP data in file ' Datastr.Info.Trial '. Skipping.']);
    return
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

if isfield(Datastr,'Force')
    
    % Do some data manipulation for the extra load
    if ~isempty(extraXLD) % Convert Nm to N and flip sign
        Datastr.Other.OtherData(:,5) = - Datastr.Other.OtherData(:,5) ./ ( 0.3*cos(Datastr.Other.OtherData(:,7)) );
        Datastr.Other.OtherData(:,6) = - Datastr.Other.OtherData(:,6) ./ ( 0.3*cos(Datastr.Other.OtherData(:,8)) );
    end
    
    % Write additional point data to the file, stored in Info field
    % Only works if you supply 'OSPelvisOffset' in extraXLD
    pelvOff = Datastr.Info.subjospelvisoffset; 
    Datastr.Other.OtherData(:,end+1) = repmat( pelvOff , [size(Datastr.Other.OtherData,1) 1]);
    Datastr.Other.OtherDataLabel{end+1} = 'subjospelvisoffset';
    
    % Create .mot file (kinetic data)
    getMot(Datastr,[subjroot '\' osFolder '\' subjroot(bckslsh(end)+1:end) trial],permVec,extraXLD,extraXLDdim);
    
else
    warning(['No field Force in trial' trial '. Skipping writing to .mot file. Unable to do ID.']);
    return;
end

% Set empty return value (to prevent saving by UI)
Datastr = [];

end