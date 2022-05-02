function [Datastr] = M1_importC3D(Datastr,markfolder,mTrans)
% gBMPDynUI markfolder=1; mTrans=1;
% 
% INPUT)
% - markfolder: string, specifying the folder containing the measurement
% trials with the marker data, relative to the subject root folder.
% Example:
% 'VZ'
% 
% -mTrans: matrix, 3x3 transformation matrix to rotate the data or swap
% dimensions. Can be empty.
% 
% OUTPUT)
% - Datastr, structure, with added fields:
% .Marker.MarkerData
% .Marker.MarkerDataLabel
% .Marker.MarkerFrameRate
% .Marker.MarkerTransform
% .Analog.AnalogData
% .Analog.AnalogFrameRate
% 
% NOTES)
% 

% Mark Vlutters - University of Twente - Jan 2016

%% Init

if isempty(mTrans)
    mTrans = eye(3);
end

%% Find file

subjrootfolder = Datastr.Info.SubjRoot;
subjtrial = Datastr.Info.Trial;

% Get all .C3D filenames in marker folder
c3dfiles = dir([subjrootfolder '\' markfolder '\*.c3d']);
c3dfilenames = cell(1,length(c3dfiles));
c3dfilenamescmp = cell(1,length(c3dfiles));
for ifile = 1:length(c3dfiles)
    c3dfilenames{ifile} = c3dfiles(ifile).name;

    % Reverse string and remove '.mat', for comparison with strncmpi
    c3dfilenamescmp{ifile} = c3dfiles(ifile).name(end-4:-1:1); 
end
clear c3dfiles;

% Check if storage filenames specified by Trials field exist. 
% If not, skip the thing
match = strncmpi( subjtrial(end:-1:1),c3dfilenamescmp,length(subjtrial) );
if sum(match) == 1
    subjmarkfile = c3dfilenames{match};
else
    warning(['Skipping C3D file ' subjtrial ': not found or multiple found.']);
    return;
end

%% Get C3D data

% Get root folder and filename
rootfolder = Datastr.Info.SubjRoot;
fullFileName = [rootfolder '\' markfolder '\' subjmarkfile];  

% Load data file
[C3DRaw] = readC3D(fullFileName);

% Do basic marker processing
C3DProc = procC3D(C3DRaw,mTrans);

% Store output
Datastr.Marker = C3DProc.Marker;
Datastr.Marker.MarkerDataTransform = mTrans;
Datastr.Analog = C3DProc.Analog;


end