function [Datastr] = A4_syncAll(Datastr,sChanAnaM,sChanAnaX,sChanTriggX,doResample, syncIMU)
% gBMPDynUI sChanAnaM=1; sChanAnaX=1; sChanTriggX=1; doResample=1; syncIMU=1;
%
%
% sChanAnaM: scalar, specifying the sync signal in the AnalogData field (marker)
% sChanAnaX: scalar, specifying the sync signal in the OtherData field (xpc-forceplates&EMG)
% sChanTriggX: scalar, specifying the sync signal in the OtherData field
% that contains the the trigger signal that was sent to the IMU system
% doResmaple: boolean, if true, attempt correction for unequal sample
% frequencies between NI box and XPC target

% sChanAnaM=1;
% sChanAnaX=1;
% sChanTriggX=2;
% doResample=false; %VZ analog is 10*100Hz(factor*markerframerate) and XPC is 1000Hz

%% Do some checks

if ~isfield(Datastr,'Other') %if no XPC sync data
    warning(['No field Other file ' Datastr.Info.Trial '. Skipping.']);
    return;
end

if ~isfield(Datastr,'Analog') % if no VZ data
    warning(['No VZ analog channel in file ' Datastr.Info.Trial '. Skipping.']);
    return
end

% Trial with analogs channels not as expected, check!
if size(Datastr.Analog.AnalogData,1)==1
    warning(['VZ analog channel in file ' Datastr.Info.Trial ' is in a row vector and not in column vectors, check (if it applies). Skipping.']);
    return
end

if ~isfield(Datastr,'IMU') % if no IMU data, synch only of XPC and VZ
    warning(['No IMU file ' Datastr.Info.Trial '. Synch only between XPC and VZ data.']);
end


%% Sync data
% Get data for syncing (here it is assumed both have the same sample freq)
sync_xpc = Datastr.Other.OtherData(:,sChanAnaX);
sync_xpc_imu = Datastr.Other.OtherData(:,sChanTriggX);
sync_c3d = Datastr.Analog.AnalogData(:,sChanAnaM);

% OPTIONAL resampling
if doResample && ~isfield(Datastr.Marker,'wasRe')
    
    % Resample marker data to XPC data
    Datastr = getResampledMarker_mod(Datastr,sync_c3d,sync_xpc);
    
    % Re-obtain data for syncing
    sync_c3d = Datastr.Analog.AnalogData(:,sChanAnaM);
end


% Get relevant sample frequencies
fs_xpc = Datastr.Other.OtherFrameRate; % >> Assumption !!!
fs_mrk = Datastr.Marker.MarkerFrameRate;
fs_imu = Datastr.IMU.IMUFrameRate;
fs_rate = fs_xpc ./ fs_mrk;
fs_rate_imu= fs_xpc ./ fs_imu;

if isfield(Datastr,'IMU') && syncIMU %Synchonize with IMU as well (IMU always within XPC length)
    flagIMU=1;
    delayIMU = find(Datastr.Other.OtherData(:,sChanTriggX),1,'first'); %Start trigger signal for IMU
    delayIMU_end = find(Datastr.Other.OtherData(:,sChanTriggX),1,'last'); %Stop trigger signal for IMU
    lengthIMU = delayIMU_end-delayIMU + 1;
else
    flagIMU=0;
end

% Find data delay (with respect to XPC!) and sync the data
datdelay = finddelay(sync_xpc,sync_c3d);



if datdelay > 0 % XPC started after VZ
    %idxStart
    idxStartXPC = 1;
    idxStartMRK = round(datdelay ./ fs_rate); %Should be corrected with 1 sample? https://nl.mathworks.com/help/signal/ref/finddelay.html
    if flagIMU %When synch with IMU, IMU always the first iDx and started always after XPC
        idxStartXPC = delayIMU;
        idxStartMRK = round((datdelay+delayIMU) ./ fs_rate); %round(((datdelay+1)+(delayIMU-1)) ./ fs_rate);
        idxStartIMU = 1;
    end
    if idxStartMRK == 0
        idxStartMRK = 1;
    end
    %idxEnd
    if find( abs(sync_c3d)>0.1 , 1,'last') == length(sync_c3d) % VZ stopped before XPC
        idxEndMRK = find(abs(sync_c3d)>0.1 ,1 ,'last') ./ fs_rate;
        idxEndXPC = find(abs(sync_c3d)>0.1 ,1 ,'last') - datdelay;
        if flagIMU %When synch with IMU
            if idxEndXPC<=delayIMU_end %stop IMU when VZ stopped
                idxEndIMU = round((idxEndXPC-delayIMU)./fs_rate_imu);
            else %stop XPC and VZ when IMU stopped
                idxEndXPC = delayIMU_end;
                idxEndIMU = length(Datastr.IMU.IMUData);
                idxEndMRK = round((delayIMU_end+datdelay)./ fs_rate) ;
            end
            %             if (length(sync_c3d)-(datdelay+delayIMU))<(lengthIMU) % % VZ stopped before XPC and IMU
            %                 idxEndMRK = round(find(abs(sync_c3d)>0.1 ,1 ,'last') ./ fs_rate); %idx of last c3d value
            %                 idxEndXPC = find(abs(sync_c3d)>0.1 ,1 ,'last') - datdelay; %idx of last c3d value (considering the delay)
            %                 idxEndIMU = round((find(abs(sync_c3d)>0.1 ,1 ,'last') - datdelay - delayIMU)./fs_rate_imu); %idx of last c3d value (considering the delay and the delay of IMU in XPC)
            %             else %(length(sync_c3d)-(datdelay+idxStartIMU))>=(idxEndIMU-idxStartIMU) % VZ stopped before XPC and after IMU
            %                 idxEndMRK = idxStartMRK + round(lengthIMU ./ fs_rate) ; %idx of StartMRK plus the length of the MIU trigg signal in 100Hz
            %                 idxEndXPC = delayIMU + lengthIMU; %idx of last IMU value
            %                 idxEndIMU = round(lengthIMU./fs_rate_imu);
            %             end
        end
        
    else % XPC stopped before VZ
        idxEndXPC = find(abs(sync_xpc)>0.1 ,1 ,'last');
        idxEndMRK = idxStartMRK + round( find(abs(sync_xpc)>0.1 ,1 ,'last') ./ fs_rate );
        %idxEndMRK = idxStartMRK + round( find(abs(sync_xpc)>0.1 ,1 ,'last') ./ fs_rate );
        
        if flagIMU %When synch with IMU   %XPC and MRK finish when IMU does
            idxEndXPC = delayIMU_end;
            idxEndIMU = length(Datastr.IMU.IMUData);
            idxEndMRK = round((delayIMU_end+datdelay)./ fs_rate) ;
        end
        
    end
    
elseif datdelay < 0 % VZ started after XPC
    %idxStart
    idxStartXPC = -datdelay+1;
    idxStartMRK = 1;
    
    if flagIMU %When synch with IMU, IMU always the first iDx and started always after XPC
        if idxStartXPC>=delayIMU % VZ started after IMU
            idxStartIMU =  round(((-datdelay+1) - delayIMU)./fs_rate_imu);
        else %VZ and XPC start when IMU does
            idxStartXPC = (-datdelay+1)+delayIMU; %-datdelay+1+delayIMU-1;
            idxStartMRK = round(delayIMU./ fs_rate);
            idxStartIMU = 1;  %Start with MRK first sample (considering delay of XPC and IMU)
            
        end
    end
    %idxEnd
    if find( abs(sync_c3d)>0.1 , 1,'last') == length(sync_c3d) % VZ stopped before XPC
        idxEndMRK = round(find(abs(sync_c3d)>0.1 ,1 ,'last') ./ fs_rate);
        idxEndXPC = find(abs(sync_c3d)>0.1 ,1 ,'last') - datdelay;
        if flagIMU %When synch with IMU
            if idxEndXPC<=delayIMU_end %stop IMU when VZ stopped
                idxEndIMU = round((idxEndXPC-delayIMU)./fs_rate_imu);
            else %stop XPC and VZ when IMU stopped
                idxEndXPC = delayIMU_end;
                idxEndIMU = length(Datastr.IMU.IMUData);
                idxEndMRK = round((delayIMU_end+datdelay)./ fs_rate) ;
            end
            %             if (length(sync_c3d)-(delayIMU-(-datdelay)))<(lengthIMU) % % VZ stopped before XPC and IMU
            %                 idxEndMRK = round(find(abs(sync_c3d)>0.1 ,1 ,'last') ./ fs_rate); %idx of last c3d value
            %                 idxEndXPC = find(abs(sync_c3d)>0.1 ,1 ,'last') - datdelay; %idx of last c3d value (considering the delay)
            %                 idxEndIMU =  round((find(abs(sync_c3d)>0.1 ,1 ,'last') - datdelay - delayIMU)./fs_rate_imu); %idx of last c3d value (considering the delay and the delay of IMU in XPC)
            %             else %(length(sync_c3d)-(datdelay+idxStartIMU))>=(idxEndIMU-idxStartIMU) % VZ stopped before XPC and after IMU
            %                 idxEndMRK = idxStartMRK + round(lengthIMU ./ fs_rate) ; %idx of StartMRK plus the length of the MIU trigg signal in 100Hz
            %                 idxEndXPC = delayIMU + lengthIMU; %idx of last IMU value
            %                 idxEndIMU =  round(lengthIMU./fs_rate_imu);
            %             end
        end
        
    else % XPC stopped before VZ
        idxEndXPC = find(abs(sync_xpc)>0.1 ,1 ,'last');
        idxEndMRK = round( (datdelay + find(abs(sync_xpc)>0.1 ,1 ,'last')) ./ fs_rate );
        if flagIMU %When synch with IMU   %XPC and MRK finish when IMU does
            idxEndXPC = delayIMU_end;
            idxEndIMU = length(Datastr.IMU.IMUData);
            idxEndMRK = round((delayIMU_end+datdelay)./ fs_rate) ;
        end
        
    end
    
end
%fprintf('#Trial %s \n',Datastr.Info.Trial)

% Store sync indices
if isfield(Datastr,'Marker')
    Datastr.Marker.MarkerSyncIdx = [idxStartMRK idxEndMRK];
    Datastr.Marker = orderfields(Datastr.Marker);
   % fprintf('Idx for IMU based on MRK: %6.f \n',(idxEndMRK-idxStartMRK ).*fs_rate./fs_rate_imu)
end
if isfield(Datastr,'Force')
    Datastr.Force.ForceSyncIdx = [idxStartXPC idxEndXPC];
    Datastr.Force = orderfields(Datastr.Force);
  %  fprintf('Idx for IMU based on XPC: %6.f \n',(idxEndXPC-idxStartXPC )./fs_rate_imu)
end


if isfield(Datastr,'EMG')
    Datastr.EMG.EMGSyncIdx = [idxStartXPC idxEndXPC];
    Datastr.EMG = orderfields(Datastr.EMG);

end
if isfield(Datastr,'Other')
    Datastr.Other.OtherSyncIdx = [idxStartXPC idxEndXPC];
    Datastr.Other = orderfields(Datastr.Other);
end

if isfield(Datastr,'IMU') && flagIMU
    Datastr.IMU.IMUSyncIdx = [idxStartIMU idxEndIMU];
    Datastr.IMU = orderfields(Datastr.IMU);
  %  fprintf('Idx for IMU based on IMU: %6.f \n',idxEndIMU-idxStartIMU)
end
%trial


