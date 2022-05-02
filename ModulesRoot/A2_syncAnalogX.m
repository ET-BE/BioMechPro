function [Datastr] = A2_syncAnalogX(Datastr,sChanAna,sChanAnaX,doResample)
% gBMPDynUI sChanAna=1; sChanAnaX=1; doResample=1;
% 
% sChanAna: scalar, specifying the sync signal in the AnalogData field
% sChanAnaX: scalar, specifying the sync signal in the OtherData field
% doResmaple: boolean, if true, attempt correction for unequal sample
% frequencies between NI box and XPC target


%% Do some checks

if ~isfield(Datastr,'Other')
    warning('No field Other found. Skipping.');
    return;
end

%% Sync data

% Get data for syncing (here it is assumed both have the same sample freq)
sync_xpc = Datastr.Other.OtherData(:,sChanAnaX);
sync_c3d = Datastr.Analog.AnalogData(:,sChanAna);


% OPTIONAL resampling
if doResample && ~isfield(Datastr.Marker,'wasRe')

    % Resample marker data to XPC data
    Datastr = getResampledMarker(Datastr,sync_c3d,sync_xpc);

    % Re-obtain data for syncing 
    sync_c3d = Datastr.Analog.AnalogData(:,sChanAna);
end


% Get relevant sample frequencies
fs_xpc = Datastr.Other.OtherFrameRate; % >> Assumption !!!
fs_mrk = Datastr.Marker.MarkerFrameRate;
fs_rate = fs_xpc ./ fs_mrk;

% Find data delay (with respect to XPC!) and sync the data
datdelay = finddelay(sync_xpc,sync_c3d);
if datdelay > 0 % XPC started after VZ
    idxStartXPC = 1;
    idxStartMRK = round(datdelay ./ fs_rate);
    if idxStartMRK == 0
        idxStartMRK = 1;
    end

    if find( abs(sync_c3d)>0.1 , 1,'last') == length(sync_c3d) % VZ stopped before XPC
        idxEndMRK = find(abs(sync_c3d)>0.1 ,1 ,'last') ./ fs_rate;
        idxEndXPC = find(abs(sync_c3d)>0.1 ,1 ,'last') - datdelay;
    else % XPC stopped before VZ
        idxEndXPC = find(abs(sync_xpc)>0.1 ,1 ,'last');
        idxEndMRK = idxStartMRK + round( find(abs(sync_xpc)>0.1 ,1 ,'last') ./ fs_rate );
    end

elseif datdelay < 0 % VZ started after XPC

    idxStartXPC = -datdelay;
    idxStartMRK = 1;

    if find( abs(sync_c3d)>0.1 , 1,'last') == length(sync_c3d) % VZ stopped before XPC
        idxEndMRK = find(abs(sync_c3d)>0.1 ,1 ,'last') ./ fs_rate;
        idxEndXPC = find(abs(sync_c3d)>0.1 ,1 ,'last') - datdelay;
    else % XPC stopped before VZ
        idxEndXPC = find(abs(sync_xpc)>0.1 ,1 ,'last');
        idxEndMRK = round( (datdelay + find(abs(sync_xpc)>0.1 ,1 ,'last')) ./ fs_rate );
    end

end

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