function [Datastr] = S21_exportEMG(Datastr, osFolder, EMGOutput, EMGFormat, norType)
% gBMPDynUI osFolder=1; EMGs_Output=1; EMGFormat=1; norType=1;
% Detrend, rectify, and low pass filter EMG
%
% INPUT)
% Datastr: structure, with at least the fields
%     .EMG.EMGData
%     .EMG.EMGFrameRate
%  EMGs_Output={'TAL','GML','RFL'}; or any EMG channel of interest
%  EMGFormat=['.mot']; %availableFileFormats=['.txt', ' .sto', ' .mot'];
%  option_fs; %=0 fs equalt to Datastr.Marker.MarkerFrameRate; =1 fs equal to selected FrameRate (defined in the code!) (default 100Hz)

%
% OUTPUT)
% No direct output
% A .mot file is created in the destination provided by filename
%
% NOTES)
% Function based on MOtion data elaboration TOolbox for
% NeuroMusculoSkeletal applications (MOtoNMS).
% Copyright (C) 2012-2014 Alice Mantoan, Monica Reggiani
%
% This function ports EMG data to a .mot file.
% With option_fs=0,  the number of
% samples written corresponds with the number of marker samples for which
% EMG data is available. It is assumed that the EMG data has AT LEAST the sample frequency of
% that of the marker data. If the EMG data has a higher sample frequency,
% some samples will be discarded.
% With option_fs=1, the data is exported at the frame rate of 100Hz by
% default (modify function accordingly for other frame rates)


%% Check input

if ~isfield(Datastr.Resample.Sych,'EMG') %if no EMG data
    warning(['No field EMG file ' Datastr.Info.Trial '. Skipping.']);
    return;
end

if isempty(EMGOutput)
    EMGOutput = Datastr.EMG.DataLabel;
end

%% Get info
% get subject folder
subjroot = Datastr.Info.SubjRoot;
trial = Datastr.Info.Trial;

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

trialOutputPath = [subjroot '\' osFolder '\' subjroot(bckslsh(end)+1:end) trial 'EMG'];

% Prepare the EMG for output
% EMG is exported at the resampling frame rate
markerFrameRate = Datastr.Resample.FrameRate;



UnNormEMG = Datastr.Resample.Sych.EMG(:, strcmp(EMGOutput, Datastr.EMG.DataLabel));

if strcmp(norType, 'dynMVC')
    dynMVC = importdata([subjroot '/dynMVCvalue.mat']);
    NormEMG = UnNormEMG./dynMVC(strcmp(EMGOutput, Datastr.EMG.DataLabel));
    Datastr.Resample.Sych.EMGExpdynNorFlag = 1;
    Datastr.Resample.Sych.EMGExpNorFlag = 0;
    
elseif strcmp(norType, 'MVC')
    MVC = importdata([subjroot '/MVCvalue.mat']);
    NormEMG = UnNormEMG./MVC(strcmp(EMGOutput, Datastr.EMG.DataLabel));
    Datastr.Resample.Sych.EMGExpdynNorFlag = 0;
    Datastr.Resample.Sych.EMGExpNorFlag = 1;
else
    NormEMG = UnNormEMG;
    Datastr.Resample.Sych.EMGExpdynNorFlag = 0;
    Datastr.Resample.Sych.EMGExpNorFlag = 0;
end

EMGtime = (0:size(NormEMG ,1)-1)'./markerFrameRate;

% ------------------------------------------------------------------------
%                            PRINT emg.txt
%--------------------------------------------------------------------------
%availableFileFormats=['.txt', ' .sto', ' .mot'];

switch EMGFormat
    
    case '.txt'
        
        printEMGtxt(trialOutputPath, EMGtime, NormEMG, EMGOutput);
        
    case {'.sto','.mot'}
        
        % fprintf('Saving filt EMG in selected format: %4.2f %\n',k/length(trialsList)*100);
        printEMGmot(trialOutputPath, EMGtime, NormEMG, EMGOutput, EMGFormat);
        
        %case ...
        %you can add here other file formats
        
    otherwise
        error('ErrorTests:convertTest', ...
            ['----------------------------------------------------------------\nWARNING: EMG Output File Format not Available!\nChoose among: [' availableFileFormats '].'])
end

cd(Datastr.Info.SubjRoot) %Make sure it is in subject folder
% Datastr=[];
end