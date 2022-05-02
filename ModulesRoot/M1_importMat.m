function [Datastr] = M1_importMat(Datastr,markfolder,mTrans)
% gBMPDynUI markfolder=1; mTrans=1;
% 
% Import .mat file into the structure, where the .mat file was obtained
% directly from the VZproject using getVZProjData.m

% - markfolder: string, specifying the folder containing the measurement
% trials with the marker data, relative to the subject root folder.
% Example:
% 'VZ'
% 
% -mTrans: matrix, 3x3 transformation matrix to rotate the data or swap
% dimensions. Can be empty.


%% Init

if isempty(mTrans)
    mTrans = eye(3);
end


%% Find file

subjrootfolder = Datastr.Info.SubjRoot;
subjtrial = [Datastr.Info.Trial 'raw']; % The 'raw' is added by getVZProjData

% Get all .mat filenames in marker folder
matfiles = dir([subjrootfolder '\' markfolder '\*.mat']);
matfilenames = cell(1,length(matfiles));
matfilenamescmp = cell(1,length(matfiles));
for ifile = 1:length(matfiles)
    matfilenames{ifile} = matfiles(ifile).name;

    % Reverse string and remove '.mat', for comparison with strncmpi
    matfilenamescmp{ifile} = matfiles(ifile).name(end-4:-1:1); 
end
clear matfiles;

% Check if storage filenames specified by Trials field exist. 
% If not, skip the thing
match = strncmpi( subjtrial(end:-1:1),matfilenamescmp,length(subjtrial) );
if sum(match) == 1
    subjmarkfile = matfilenames{match};
else
    try
        warning(['Skipping Mat file ' subjtrial ': not found or multiple found.']);
        return;
    catch
        error('A1_importMat:nofile',['Skipping Mat file ' subjtrial ': not found.']);
    end
end

%% Get .mat data

% Get root folder and filename
rootfolder = Datastr.Info.SubjRoot;
fullFileName = [rootfolder '\' markfolder '\' subjmarkfile];  

% Load data file
load(fullFileName);

% Remove jumps (disabled by default)
% C3DRaw = removeMarkerJumps(C3DRaw);

% Do basic marker processing
C3DProc = procC3D(C3DRaw,mTrans);

% Store output
Datastr.Marker = C3DProc.Marker;
Datastr.Marker.MarkerDataTransform = mTrans;
Datastr.Analog = C3DProc.Analog;


end