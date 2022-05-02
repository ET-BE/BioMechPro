function [Datastr] = S3_importINSOLE(Datastr,insoleFolder,leftID,rightID, frameRate)
% gBMPDynUI insoleFolder=1; leftID=1; rightID=1; frameRate=1;
% 
% Load the insole text file, process, and sve them into the Datastr
% 
% INPUT)
% Datastr: structure, loaded marker/IMU data
%
% insoleFolder: string, specifying the folder containing the measurement
% trials with the insole data, relative to the subject root folder.
% Example:
% 'Insole'
%
% leftID: integer, the number ID of the left insole data. Default: 0
% rightID: integer, ther number ID of the right insole data. Default: 1
% frameRate: integer, the sampling rate of insole recording.
% 
% OUTPUT)
% Datastr: structure with no new added fields

%% Init

if isempty(leftID)
    leftID = 0;
end

if isempty(rightID)
    rightID = 1;
end

if isempty(frameRate)
    frameRate = 100;
end

%% Find file

subjrootfolder = Datastr.Info.SubjRoot;
subjtrial = Datastr.Info.Trial;

% Get all .txt filenames in insole folder
insolefiles = dir([subjrootfolder '\' insoleFolder '\*.txt']);
insolefilenames = cell(1,length(insolefiles));
insolefilenamescmp = cell(1,length(insolefiles));
for ifile = 1:length(insolefiles)
    insolefilenames{ifile} = insolefiles(ifile).name;

    % Reverse string and remove '.mat', for comparison with strncmpi
    insolefilenamescmp{ifile} = insolefiles(ifile).name(end-4:-1:1); 
end
clear insolefiles;

% Check if storage filenames specified by Trials field exist. 
% If not, skip the thing
match = strncmpi( subjtrial(end:-1:1), insolefilenamescmp, length(subjtrial));
if sum(match) == 1
    subjinsolefile = insolefilenames{match};
else
    warning(['Skipping insole file ' subjtrial ': not found or multiple found.']);
    return;
end

%% Get Insole data

% get root folder and filename
rootfolder = Datastr.Info.SubjRoot;
fullFileName = [rootfolder '\' insoleFolder '\' subjinsolefile];

% load data file
[insoleRaw] = importdata(fullFileName);
[insoleProc] = procInsole(insoleRaw, leftID, rightID);

% save data
Datastr.Insole.FrameRate = frameRate;
Datastr.Insole.DataLabel = string(insoleRaw.textdata);
Datastr.Insole.DataLabel(2) = [];

Datastr.Insole.Left = insoleProc.left;
Datastr.Insole.Right = insoleProc.right;
end


