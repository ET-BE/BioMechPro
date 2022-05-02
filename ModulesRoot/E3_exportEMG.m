function [Datastr] = E3_exportEMG(Datastr,EMGs_Output,EMGFormat,option_fs)
% gBMPDynUI EMGs_Output=1; EMGFormat=1; option_fs=1;
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

if ~isfield(Datastr,'EMG') %if no EMG data
    warning(['No field EMG file ' Datastr.Info.Trial '. Skipping.']);
    return;
end

% Check if sync indices exist for MarkerData
if isfield(Datastr,'Marker')
    if ~isfield(Datastr.Marker,'MarkerSyncIdx')
        warning('Marker:syncIdx',['No sync indices found in ' Datastr.Info.Trial '. Assuming synchronized Marker data']);
        markerSyncIdx = [1 size(Datastr.Marker.MarkerData,1)];
    else
        markerSyncIdx = Datastr.Marker.MarkerSyncIdx;
    end
end

% Check if sync indices exist for IMUData
if isfield(Datastr,'IMU')
    if ~isfield(Datastr.IMU,'IMUSyncIdx')
        warning('EMG:syncIdx',['No sync indices found in ' Datastr.Info.Trial '. Assuming synchronized IMU data']);
        IMUSyncIdx = [1 size(Datastr.IMU.IMUData,1)];
    else
        IMUSyncIdx = Datastr.IMU.IMUSyncIdx;
    end
end

% Check if sync indices exist for EMGData
if ~isfield(Datastr.EMG,'EMGSyncIdx')
    warning('EMG:syncIdx',['No sync indices found in ' Datastr.Info.Trial '. Assuming synchronized EMG data']);
    
    EMGSyncIdx = [1 size(Datastr.EMG.EMGData,1)];
else
    EMGSyncIdx = Datastr.EMG.EMGSyncIdx;
end

%% Do Stuff
%Input BioMechPro

subjroot=Datastr.Info.SubjRoot;
foo = strfind(subjroot,'\');
infilename = subjroot(foo(end)+1:end);

trialOutputPathcell= mkOutputDir([subjroot filesep 'EMG'],{[infilename Datastr.Info.Trial]});
trialOutputPath=trialOutputPathcell{1};

%Depending on option_fs, EMG is exported in different frame rate

switch option_fs
    case 0
        nMarkSmpl = markerSyncIdx(end) - markerSyncIdx(1) + 1;
        EMGIdx = round(linspace(EMGSyncIdx(1),EMGSyncIdx(end),nMarkSmpl));

        % Frame rates
        markerFrameRate = Datastr.Marker.MarkerFrameRate;

        for icol=1:length(EMGs_Output)
            k=find(strcmp(EMGs_Output{icol},Datastr.EMG.EMGDataLabel));
            NormEMG(:,icol)= Datastr.EMG.normEMGData(EMGIdx,k);
        end
        EMGtime = (0:size(NormEMG ,1)-1)'./markerFrameRate;
    case 1
        % Frame rates
        EMGFrameRate = Datastr.EMG.EMGFrameRate;
        StandardFrameRate=100; % Modify to frame rate desired
        [N,D]=rat(StandardFrameRate/EMGFrameRate);
        try
            normEMGData = resample(Datastr.EMG.normEMGData,N,D);
        catch
            disp('error')
        end
        for icol=1:length(EMGs_Output)
            k=find(strcmp(EMGs_Output{icol},Datastr.EMG.EMGDataLabel));
            NormEMG(:,icol)= normEMGData(:,k);
        end
        EMGtime = (0:size(NormEMG ,1)-1)'./StandardFrameRate;%
end


% ------------------------------------------------------------------------
%                            PRINT emg.txt
%--------------------------------------------------------------------------
%availableFileFormats=['.txt', ' .sto', ' .mot'];

switch EMGFormat
    
    case '.txt'
        
        printEMGtxt(trialOutputPath,EMGtime,NormEMG,EMGs_Output);
        
    case {'.sto','.mot'}
        
        % fprintf('Saving filt EMG in selected format: %4.2f %\n',k/length(trialsList)*100);
        printEMGmot(trialOutputPath,EMGtime,NormEMG,EMGs_Output, EMGFormat);
        
        %case ...
        %you can add here other file formats
        
    otherwise
        error('ErrorTests:convertTest', ...
            ['----------------------------------------------------------------\nWARNING: EMG Output File Format not Available!\nChoose among: [' availableFileFormats '].'])
end

cd(Datastr.Info.SubjRoot) %Make sure it is in subject folder
Datastr=[];
end