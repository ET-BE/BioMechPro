function [Datastr] = A7_CropEMGIMU(Datastr,hasIMU)
% gBMPDynUI hasIMU=1;

%% Do some checks

if isfield(Datastr,'Window') %if no XPC sync data
    if isfield(Datastr.Window,'wasCrop')
        warning(['Data already cropped in ' Datastr.Info.Trial '.']);
        return
    end
end

if ~isfield(Datastr,'IMU') && ~hasIMU% if no IMU data
    warning(['No IMU file ' Datastr.Info.Trial '. Skipping.']);
    return
end

if ~isfield(Datastr,'EMG') % if no EMG data
    warning(['No EMG file ' Datastr.Info.Trial '. Skipping.']);
    return
end

if ~isfield(Datastr.EMG,'EMGSyncIdx') ||  ~isfield(Datastr.IMU,'IMUSyncIdx')
    warning(['No sync indices found in' infilename '. Assuming synchronized EMG-IMU data']);
    
    EMGSyncIdx = [1 size(Datastr.EMG.EMGData,1)];
    imuSyncIdx = [1 size(Datastr.IMU.IMUData,1)];
else
    
    EMGSyncIdx = Datastr.EMG.EMGSyncIdx;
    imuSyncIdx = Datastr.IMU.IMUSyncIdx;
end

%% 
subjroot=Datastr.Info.SubjRoot;
itrial=Datastr.Info.Trial;
%Get all signals in the synchronized window
otherdata=Datastr.Other.OtherData(EMGSyncIdx(1):EMGSyncIdx(end),:);
emgdata=Datastr.EMG.EMGData(EMGSyncIdx(1):EMGSyncIdx(end),:);
normemgdata=Datastr.EMG.normEMGData(EMGSyncIdx(1):EMGSyncIdx(end),:);
filtemgdata=Datastr.EMG.filtEMGData(EMGSyncIdx(1):EMGSyncIdx(end),:);
imudata=Datastr.IMU.IMUData(imuSyncIdx(1):imuSyncIdx(end),:);
contactimudata=Datastr.IMU.IMUContactRFoot(imuSyncIdx(1):imuSyncIdx(end),:);

fs_xpc = Datastr.Other.OtherFrameRate;
fs_imu = Datastr.IMU.IMUFrameRate;


%Plotting
figure('name', ['Trial: ' itrial] ,...
    'color',[1 1 1],...
    'units', 'normalized',...
    'outerposition',[0 0 1 1]);
% 1) JointAngle synch
%             plot(idxStartM*Ts,ones(1,length(idxStartM))*700,'r*',idxStopM*Ts,ones(1,length(idxStopM))*700,'o')
%             text(idxStartM*Ts-.5,ones(1,length(idxStartM))*700,num2str((1:numel(idxStartM))'));
imuData = Datastr.IMU.IMUData(:,:);
t=(0:size(imuData,1)-1)'./fs_imu;
h=plot(t,imuData(:,[10,13])); %kne and ankle angles
xlabel('time(s)')
ylabel(['Joint angle (' char(176) ')'])
% title(['Trial ' num2str(itrial)])
legend(h,[{'Knee'} {'Ankle'}],'Location','southwest');

% Time to cut
if isfield(Datastr,'Window')
    %Get sample at which starts
    tstart=Datastr.Window.TimeinMrk(1);
    %Get sample at which ends
    tend=Datastr.Window.TimeinMrk(2);
    idxStartO = round( tstart * fs_xpc )  ;
    idxStopO = round( tend * fs_xpc )  ;
    idxStartI = round( tstart * fs_imu )  ;
    idxStopI = round( tend * fs_imu )  ;
    
else
    %Select times
    prompt = {'Enter start time','Enter end time'};
    title_w = 'Select window time';
    dims = [1 35];
    definput = {'10','20'};
    answer = inputdlg(prompt,title_w,dims,definput);
    %Get time at which starts
    if strcmp(answer{1},'start')
        tstart=0;
        idxStartM = 1;
        idxStartO = 1;
        idxStartI = 1 ;
    else
        tstart=str2num(answer{1});
        idxStartO = round( tstart * fs_xpc )  ;
        idxStartI = round( tstart * fs_imu )  ;
    end
    
    
    
    %Get time at which ends
    if strcmp(answer{2},'end')
        tend=size(imudata,1)/fs_imu;
        idxStopO = size(otherdata,1)  ;
        idxStopI = size(imudata,1)  ;
    else
        tend=str2num(answer{2});
        idxStopO = round( tend * fs_xpc )  ;
        idxStopI = round( tend * fs_imu )  ;
    end
    Datastr.Window.TimeinMrk=[tstart tend];
    %     save([subjroot  '_' num2str(itrial) '.mat'],'Datastr')
    
end

% From here on data that is not within the synchronized window is deleted
Datastr.Other.OtherData=otherdata(idxStartO:idxStopO,:,:);
Datastr.EMG.EMGData=emgdata(idxStartO:idxStopO,:);
Datastr.EMG.normEMGData=normemgdata(idxStartO:idxStopO,:);
Datastr.EMG.filtEMGData=filtemgdata(idxStartO:idxStopO,:);
Datastr.IMU.IMUData=imudata(idxStartI:idxStopI,:);
Datastr.IMU.IMUContactRFoot=contactimudata(idxStartI:idxStopI,:);

%Renew idx of Synch
Datastr.Other.OtherSyncIdx = [1 size(Datastr.Other.OtherData,1)];
Datastr.EMG.EMGSyncIdx = [1 size(Datastr.EMG.EMGData,1)];
Datastr.IMU.IMUSyncIdx = [1 size(Datastr.IMU.IMUData,1)];

Datastr.Window.wasCrop=1;

