function [Datastr] = S1_importMatQTM(Datastr, matQTMfolder, fSign, ...
                     fTrans,eSign, AnaEMGInd)
% gBMPDynUI matQTMfolder=1; fSign=1; fTrans=1; emgSign=1; AnaEMGInd=1; 
% , mvcTrialName, fType, fOrd, mfCut, ffCut, efCut
% mvcFolder=1; fType=1; fOrd=1; mfCut=1; ffCut=1; efCut=1; 
% 
% Import .mat file into the structure, where the .mat file was obtained
% by exporting from QTM (the unit for the marker data is mm (cannot be 
% changed (?), this cannot be automatically retrieved)

% - matQTMfolder: string, specifying the folder containing the measurement
% trials with the marker data (may contians force and emg data),
% relative to the subject root folder.
% Example:
% 'MAT'
%
% - fSign: string, specifying whether the treadmill force data are
% included or not.
%
% - fTrans: vector, change the force data order into the same as marker,
% also the same as the OpenSim model setup:
% For example, for the lops lab, the fTrans is [2, 3, 1]
%
%        ----------- -----------
%        |    L    |x|    R    |
%        |         |^|         |
%        |    |    |||         |
%        |    |    |||         |
%        |    |    |*|------> z|
%        |    V    |y|         |
%        |   belt  | |         |
%        |direction| |         |
%        |         | |         |
%        ----------- -----------
% 
% -emgSign: string, specifying whether the treadmill force data are
% included or not.
%
% - AnaEMGInd: vector, the index vector of the EMG recordings in the Analog
% channels.
%
% - mvcTrialName: string, the MVC trial for EMG normalization.
%
% - fType: string, filter type for marker, 
%
% - fOrd: integer, the order of filter
%
% - mfCut: integer, the cutoff frequency of the filter for markers;
%
% - ffCut: integer, the cutoff freqency of the filter for forces;
%
% - efCut: 1x2 vector, the cutoff freqencies of the filters for EMGs
%
% Adjusted: Huawei Wang - University of Twente - May 2021

%% Init inputs transform vector

if isempty(fTrans)
    fTrans = [1, 2, 3];
end


%% Find file
subjrootfolder = Datastr.Info.SubjRoot;
subjtrial = [Datastr.Info.Trial];

% Get all .mat filenames in marker folder
matfiles = dir([subjrootfolder '\' matQTMfolder '\*.mat']);
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
fullFileName = [rootfolder '\' matQTMfolder '\' subjmarkfile];  

% Load data file
QTMData=importdata(fullFileName);

% Put data into correct structure (Datastr that is used for further processing)
% Marker data:
DatastrRaw.Marker.FrameRate=QTMData.FrameRate;
DatastrRaw.Marker.DataLabel=QTMData.Trajectories.Labeled.Labels;
DatastrRaw.Marker.MarkerData=squeeze(QTMData.Trajectories.Labeled.Data(:,1:3,:))/1000;
DatastrRaw.Marker.units = 'm';

% Store marker data into Datastr
Datastr.Marker = DatastrRaw.Marker;

if strcmp(fSign, 'True')
    % Force data: Label
    DatastrRaw.Force.DataLabel = Datastr.Info.ForceLabel;
    
    % save data based on the order of force label and transform vector
    % ["Fx", "Fy", "Fz", "Mx", "My", "Mz",  "CoPx", "CoPy", "CoPz"]  
    for i = 1:size(QTMData.Force, 2)
        if strcmp(char(QTMData.Force(i).ForcePlateName), 'Force-plate Right')
            DatastrRaw.Force.RightForceData = [sign(fTrans).*QTMData.Force(i).Force(abs(fTrans), :)',...
                                      sign(fTrans).*QTMData.Force(i).Moment(abs(fTrans), :)',...
                                      sign(fTrans).*QTMData.Force(i).COP(abs(fTrans), :)'];
            DatastrRaw.Force.RightPlateLocation = QTMData.Force(i).ForcePlateLocation;
            DatastrRaw.Force.FrameRate = QTMData.Force(i).Frequency;
        elseif strcmp(char(QTMData.Force(i).ForcePlateName), 'Force-plate Left')
            DatastrRaw.Force.LeftForceData = [sign(fTrans).*QTMData.Force(i).Force(abs(fTrans), :)',...
                                     sign(fTrans).*QTMData.Force(i).Moment(abs(fTrans), :)',...
                                     sign(fTrans).*QTMData.Force(i).COP(abs(fTrans), :)'];
            DatastrRaw.Force.LeftPlateLocation = QTMData.Force(i).ForcePlateLocation;
        end
    end
    
    Datastr.Force.ForceDataTransform = fTrans;
    Datastr.Force = DatastrRaw.Force;
    
end

if strcmp(eSign, 'True')
    % extract EMG data
    DatastrRaw.EMG.FrameRate=QTMData.Analog.Frequency;
    DatastrRaw.EMG.DataLabel=Datastr.Info.EMGLabel;
    DatastrRaw.EMG.Channels =QTMData.Analog.Data(AnaEMGInd, :)';
    
    Datastr.EMG = DatastrRaw.EMG;
end

end