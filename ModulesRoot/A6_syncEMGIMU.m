function [Datastr] = A6_syncEMGIMU(Datastr,sChanTriggX)
% gBMPDynUI sChanTriggX=1;
% 
% sChanTriggX: scalar, specifying the sync signal in the OtherData field
% that contains the the trigger signal that was sent to the IMU system

%% Do some checks

if ~isfield(Datastr,'Other')
    warning(['No field Other file ' Datastr.Info.Trial '. Skipping.']);
    return;
end
if ~isfield(Datastr,'IMU') % if no IMU data, synch only of XPC and VZ
    warning(['No IMU file ' Datastr.Info.Trial '. Skipping.']);
end

%% Sync data

%Synch EMG according to when rising and falling edges were detected
idxStartXPC = find(Datastr.Other.OtherData(:,sChanTriggX),1,'first'); %Start of trigger signal for IMU in XPC data
idxEndXPC = find(Datastr.Other.OtherData(:,sChanTriggX),1,'last'); %Stop of trigger signal for IMU in XPC data


if isfield(Datastr,'IMU')
    Datastr.IMU.IMUSyncIdx = [1 size(Datastr.IMU.IMUData,1)];
    Datastr.IMU = orderfields(Datastr.IMU);
end
if isfield(Datastr,'EMG')
    Datastr.EMG.EMGSyncIdx = [idxStartXPC idxEndXPC];
    Datastr.EMG = orderfields(Datastr.EMG);
end
if isfield(Datastr,'Other')
    Datastr.Other.OtherSyncIdx = [idxStartXPC idxEndXPC];
    Datastr.Other = orderfields(Datastr.Other);
end

end