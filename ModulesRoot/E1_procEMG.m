function [Datastr] = E1_procEMG(Datastr,doBStop,fOrd,fCut)
% gBMPDynUI doBStop=1; fOrd=1; fCut=1;
% 
% Detrend, rectify, and low pass filter EMG
% 
% INPUT)
% Datastr: structure, with at least the fields
%     .EMG.EMGData
%     .EMG.EMGFrameRate
% 
% OUTPUT)
% Datastr: structure with no new added fields
% 
% NOTES)
% No zero phase filter (filtfilt) is used here

%% Check input

if ~isfield(Datastr,'EMG')
    return;
end

%% Do Stuff

% Get sample frequency
fs_emg = Datastr.EMG.EMGFrameRate;

% Get EMG
emgData = Datastr.EMG.EMGData;

% Detrend
emgData = emgData - repmat( mean(emgData,1) , [size(emgData,1) 1] );

if doBStop
    % Create 50 Hz bandstop filter
    [bstop,astop] = butter(1,2.*[48 52]./fs_emg,'stop');

    % Filter for 50 Hz
    emgData = filter(bstop,astop,emgData);
end

% Rectify
emgData = abs(emgData);

% Create lowpass filter
[blow,alow] = butter(fOrd,2.*fCut./fs_emg,'low');

% Filter EMG
emgData = filter(blow,alow,emgData);

% Store data
Datastr.EMG.EMGData = emgData;


end