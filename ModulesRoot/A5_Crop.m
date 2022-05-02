function [Datastr] = A5_Crop(Datastr,hasIMU,hasMarker,hasForce)
% gBMPDynUI hasIMU=1; hasMarker=1; hasForce=1;

%% Do some checks

if isfield(Datastr,'Window') %if no XPC sync data
    if isfield(Datastr.Window,'wasCrop')
        warning(['Data already cropped in ' Datastr.Info.Trial '.']);
        return
    end
end

if ~isfield(Datastr,'Marker') %if no VZ data
    warning(['No field Marker file ' Datastr.Info.Trial '. Might crash.']);
    if ~hasMarker
        return
    end
end

if ~isfield(Datastr,'Force') % if no XPC data
    warning(['No field Force file ' Datastr.Info.Trial '. Might crash.']);
    if ~hasForce
        return
    end
end


if ~isfield(Datastr,'IMU') % if no IMU data
    warning(['No IMU file ' Datastr.Info.Trial '. Might crash.']);
    if ~hasIMU
        return
    end
end

if ~isfield(Datastr.Marker,'MarkerSyncIdx') || ~isfield(Datastr.Force,'ForceSyncIdx') || ~isfield(Datastr.IMU,'IMUSyncIdx')
    warning('getMot:syncIdx',['No sync indices found in' infilename '. Assuming synchronized marker-force-IMU data']);
    
    markerSyncIdx = [1 size(Datastr.Marker.MarkerData,1)];
    forceSyncIdx = [1 size(Datastr.Force.ForceData,1)];
    otherSyncIdx = [1 size(Datastr.Other.OtherData,1)];
    imuSyncIdx = [1 size(Datastr.IMU.IMUData,1)];
else
    markerSyncIdx = Datastr.Marker.MarkerSyncIdx;
    forceSyncIdx = Datastr.Force.ForceSyncIdx;
    otherSyncIdx = Datastr.Other.OtherSyncIdx;
    imuSyncIdx = Datastr.IMU.IMUSyncIdx;
end

%%
subjroot=Datastr.Info.SubjRoot;
itrial=Datastr.Info.Trial;


markerdata=Datastr.Marker.MarkerData(markerSyncIdx(1):markerSyncIdx(end),:,:);
try
probeddata=Datastr.Marker.ProbedData(markerSyncIdx(1):markerSyncIdx(end),:,:);
catch
end

%Get all signals in the synchronized window
forcedata=Datastr.Force.ForceData(forceSyncIdx(1):forceSyncIdx(end),:,:);
otherdata=Datastr.Other.OtherData(otherSyncIdx(1):otherSyncIdx(end),:);
emgdata=Datastr.EMG.EMGData(forceSyncIdx(1):forceSyncIdx(end),:);
normemgdata=Datastr.EMG.normEMGData(forceSyncIdx(1):forceSyncIdx(end),:);
filtemgdata=Datastr.EMG.filtEMGData(forceSyncIdx(1):forceSyncIdx(end),:);
imudata=Datastr.IMU.IMUData(imuSyncIdx(1):imuSyncIdx(end),:);
imucontactdata=Datastr.IMU.IMUContactRFoot(imuSyncIdx(1):imuSyncIdx(end),:);

fs_mark = Datastr.Marker.MarkerFrameRate;
fs_xpc = Datastr.Other.OtherFrameRate;
fs_imu = Datastr.IMU.IMUFrameRate;


%Plotting force and joint angles to decide where to cut
f=figure(str2num(itrial));
p = uipanel('Parent',f,'BorderType','none');
p.Title = ['Trial ' itrial];
p.TitlePosition = 'centertop';
p.FontSize = 12;
p.FontWeight = 'bold';
f.Units='normalized';
f.OuterPosition=[0 0 1 1];
sub1 = subplot(5,1,1:5,'Parent',p);

% 1) Force synch
nMarkSmpl = markerSyncIdx(end) - markerSyncIdx(1) + 1;
forceIdx = round(linspace(forceSyncIdx(1),forceSyncIdx(end),nMarkSmpl));
forceData = Datastr.Force.ForceData(forceIdx,:);
t=(0:size(forceData,1)-1)'./fs_mark;
yyaxis left
ch_force=9; %Fz right
%             patch([idxStartM(i_left)*Ts idxStartM(i_right)*Ts idxStartM(i_right)*Ts idxStartM(i_left)*Ts], [min(forceData(:,ch_force))*[1 1] max(forceData(:,ch_force))*[1 1]], [0.9 0.9 0.9])
%             hold on
h_1=plot(sub1,t,forceData(:,ch_force));
ylabel('Force')
hold on
% 2)JointAngle synch
%             plot(idxStartM*Ts,ones(1,length(idxStartM))*700,'r*',idxStopM*Ts,ones(1,length(idxStopM))*700,'o')
%             text(idxStartM*Ts-.5,ones(1,length(idxStartM))*700,num2str((1:numel(idxStartM))'));
imuIdx = round(linspace(imuSyncIdx(1),imuSyncIdx(end),nMarkSmpl));
imuData = Datastr.IMU.IMUData(imuIdx,:);
t=(0:size(imuData,1)-1)'./fs_mark;
yyaxis right
h_2=plot(sub1,t,imuData(:,[10,13],1)); %knee and ankle joint angles in IMU, flexion extension
xlabel('time(s)')
ylabel(['Joint angle (' char(176) ')'])
legend(h_2,[{'Knee'} {'Ankle'}],'Location','southwest');

% Time to cut
if isfield(Datastr,'Window')
    %Get sample at which starts
    tstart=Datastr.Window.TimeinMrk(1);
    %Get sample at which ends
    tend=Datastr.Window.TimeinMrk(2);
    
    idxStartM = round( tstart * fs_mark )  ;
    idxStopM = round( tend * fs_mark )  ;
    idxStartO = round( tstart * fs_xpc )  ;
    idxStopO = round( tend * fs_xpc )  ;
    idxStartI = round( tstart * fs_imu )  ;
    idxStopI = round( tend * fs_imu )  ;
    
else
    %Select times from emerging window
    prompt = {'Enter start time','Enter end time'};
    title_w = 'Select window time';
    dims = [1 35];
    definput = {'10','20'};
    answer = inputdlg(prompt,title_w,dims,definput);
    %Get time at which starts
    if strcmp(answer{1},'start')
        tstart=0.1; %Always remove some frames at the beginning: 0.1s
        idxStartM = round( tstart * fs_mark )  ;
        idxStartO = round( tstart * fs_xpc )  ;
        idxStartI = round( tstart * fs_imu )  ;
    else
        tstart=str2num(answer{1});
        idxStartM = round( tstart * fs_mark )  ;
        idxStartO = round( tstart * fs_xpc )  ;
        idxStartI = round( tstart * fs_imu )  ;
    end   
    %Get time at which ends
    if strcmp(answer{2},'end')
        tend=(size(markerdata,1)-10)/fs_mark; %%Always remove some frames at the end: 10frames*0.01=0.1s
        idxStopM = round( tend * fs_mark )  ;
        idxStopO = round( tend * fs_xpc )  ;
        idxStopI = round( tend * fs_imu )  ;
    else
        tend=str2num(answer{2});
        idxStopM = round( tend * fs_mark )  ;
        idxStopO = round( tend * fs_xpc )  ;
        idxStopI = round( tend * fs_imu )  ;
    end
    Datastr.Window.TimeinMrk=[tstart tend];
    %     save([subjroot  '_' num2str(itrial) '.mat'],'Datastr')
    
end


% From here on data that is not within the synchronized window is deleted
Datastr.Marker.MarkerData= markerdata(idxStartM:idxStopM,:,:);
try
Datastr.Marker.ProbedData= probeddata(idxStartM:idxStopM,:,:);
catch
end
Datastr.Force.ForceData=forcedata(idxStartO:idxStopO,:,:);
Datastr.Other.OtherData=otherdata(idxStartO:idxStopO,:,:);
Datastr.EMG.EMGData=emgdata(idxStartO:idxStopO,:);
Datastr.EMG.normEMGData=normemgdata(idxStartO:idxStopO,:);
Datastr.EMG.filtEMGData=filtemgdata(idxStartO:idxStopO,:);
Datastr.IMU.IMUData=imudata(idxStartI:idxStopI,:);
Datastr.IMU.IMUContactRFoot =imucontactdata(idxStartI:idxStopI,:);

%Renew idx of Synch
Datastr.Marker.MarkerSyncIdx = [1 size(Datastr.Marker.MarkerData,1)];
Datastr.Force.ForceSyncIdx = [1 size(Datastr.Force.ForceData,1)];
Datastr.Other.OtherSyncIdx = [1 size(Datastr.Other.OtherData,1)];
Datastr.EMG.EMGSyncIdx = [1 size(Datastr.EMG.EMGData,1)];
Datastr.IMU.IMUSyncIdx = [1 size(Datastr.IMU.IMUData,1)];


Datastr.Window.wasCrop=1;

