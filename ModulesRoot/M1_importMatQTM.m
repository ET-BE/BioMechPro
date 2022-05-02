function [Datastr] = M1_importMatQTM(Datastr,markfolder,mTrans)
% gBMPDynUI markfolder=1; mTrans=1;
% 
% Import .mat file into the structure, where the .mat file was obtained
% by exporting from QTM (the unit for the marker data is mm (cannot be 
% changed (?), this cannot be automatically retrieved)

% - markfolder: string, specifying the folder containing the measurement
% trials with the marker data, relative to the subject root folder.
% Example:
% 'MAT'
% 
% -mTrans: matrix, 3x3 transformation matrix to rotate the data or swap
% dimensions. Can be empty.


%% Init

if isempty(mTrans)
    mTrans = eye(3);
end


%% Find file

subjrootfolder = Datastr.Info.SubjRoot;
subjtrial = [Datastr.Info.Trial]; % The 'raw' is added by getVZProjData

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
QTMData=importdata(fullFileName);

% Put data into correct structure (Datastr that is used for further processing)
% Marker data:
DatastrRaw.Marker.VideoFrameRate=QTMData.FrameRate;
DatastrRaw.Marker.MarkerDataLabel=QTMData.Trajectories.Labeled.Labels;
DatastrRaw.Marker.MarkerData=permute(QTMData.Trajectories.Labeled.Data(:,1:3,:),[3,1,2]);
% Analog data:
DatastrRaw.Analog.AnalogFrameRate=QTMData.Analog.Frequency;
DatastrRaw.Analog.AnalogDataLabel=QTMData.Analog.Labels;
DatastrRaw.Analog.AnalogData=QTMData.Analog.Data';
% Parameters (this is now still based on format used for C3D files and to be able to use procC3D, might be coded in a different way in the future):
DatastrRaw.ParameterGroup.units='mm';
DatastrRaw.ParameterGroup.labels=QTMData.Trajectories.Labeled.Labels;

% Remove jumps (disabled by default)
% C3DRaw = removeMarkerJumps(C3DRaw);

% Do basic marker processing
DatastrProc = procMatQTM(DatastrRaw,mTrans);

% Store output
Datastr.Marker = DatastrProc.Marker;
Datastr.Marker.MarkerDataTransform = mTrans;
Datastr.Analog = DatastrProc.Analog;


end