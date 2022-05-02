function [Datastr] = A3_syncPulse(Datastr,sChanTriggM,sChanTriggX)
% gBMPDynUI sChanTriggM=1; sChanTriggX=1; doResample=1;
% 
% Synchronization of XPC and VZ data via a rising and falling edge
% signal
% sChanTriggM: scalar, specifying the sync signal (trigger from Simulink) in the AnalogData field (marker)
% sChanTriggX: scalar, specifying the sync signal (trigger from Simulink) in the OtherData field (xpc-forceplates&EMG)


%% Do some checks

if ~isfield(Datastr,'Other')
    warning(['No field Other file ' Datastr.Info.Trial '. Skipping.']);
    return;
end
if ~isfield(Datastr,'Analog')
    warning(['No field Analog file ' Datastr.Info.Trial '. Skipping.']);
    return;
end

%% Sync data


idxStartXPC = find(Datastr.Other.OtherData(:,sChanTriggX),1,'first'); %Start of trigger signal for IMU in XPC data
idxEndXPC = find(Datastr.Other.OtherData(:,sChanTriggX),1,'last'); %Stop of trigger signal for IMU in XPC data
idxStartMRK = find(Datastr.Analog.AnalogData(:,sChanTriggM),1,'first'); %Start of trigger signal for IMU in VZ
idxEndMRK = find(Datastr.Analog.AnalogData(:,sChanTriggM),1,'last'); %Stop of trigger signal for IMU in VZ

% Store sync indices
Datastr.Marker.MarkerSyncIdx = [idxStartMRK idxEndMRK];
Datastr.Marker = orderfields(Datastr.Marker);

if isfield(Datastr,'Force')
    Datastr.Force.ForceSyncIdx = [idxStartXPC idxEndXPC];
    Datastr.Force = orderfields(Datastr.Force);
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