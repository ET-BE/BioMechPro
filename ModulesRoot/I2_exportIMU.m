% /*
%     WARNING: COPYRIGHT (C) 2018 XSENS TECHNOLOGIES OR SUBSIDIARIES
%     WORLDWIDE. ALL RIGHTS RESERVED. THIS FILE AND THE SOURCE CODE IT
%     CONTAINS (AND/OR THE BINARY CODE FILES FOUND IN THE SAME FOLDER THAT
%     CONTAINS THIS FILE) AND ALL RELATED SOFTWARE (COLLECTIVELY, "CODE")
%     ARE SUBJECT TO A RESTRICTED LICENSE AGREEMENT ("AGREEMENT") BETWEEN
%     XSENS AS LICENSOR AND THE AUTHORIZED LICENSEE UNDER THE AGREEMENT.
%     THE CODE MUST BE USED SOLELY WITH XSENS PRODUCTS INCORPORATED INTO
%     LICENSEE PRODUCTS IN ACCORDANCE WITH THE AGREEMENT. ANY USE,
%     MODIFICATION, COPYING OR DISTRIBUTION OF THE CODE IS STRICTLY
%     PROHIBITED UNLESS EXPRESSLY AUTHORIZED BY THE AGREEMENT. IF YOU ARE
%     NOT AN AUTHORIZED USER OF THE CODE IN ACCORDANCE WITH THE AGREEMENT,
%     YOU MUST STOP USING OR VIEWING THE CODE NOW, REMOVE ANY COPIES OF THE
%     CODE FROM YOUR COMPUTER AND NOTIFY XSENS IMMEDIATELY BY EMAIL TO
%     INFO@XSENS.COM. ANY COPIES OR DERIVATIVES OF THE CODE (IN WHOLE OR IN
%     PART) IN SOURCE CODE FORM THAT ARE PERMITTED BY THE AGREEMENT MUST
%     RETAIN THE ABOVE COPYRIGHT NOTICE AND THIS PARAGRAPH IN ITS ENTIRETY,
%     AS REQUIRED BY THE AGREEMENT.
% */


function [Datastr] = I2_exportIMU(Datastr,imufolder,option_fs)
% gBMPDynUI imufolder=1; option_fs=1;
% imufolder=['IMU'];
% Create a .mot file with join angles for OpenSim input 
% INPUT)
% - imufolder: string, specifying the folder containing the measurement
% trials with the IMU data, relative to the subject root folder.
% Example:
% 'IMU'
%
% option_fs; %=0 Datastr.Marker.MarkerFrameRate; =1 Selected FrameRate (default 100Hz)
%
% OUTPUT)
% No direct output
% A .mot file is created in the destination provided by filename
%
% NOTES)
%
% This function ports IMU data to a .mot file.
% With option_fs=0,  the number of 
% samples written corresponds with the number of marker samples for which 
% IMU data is available. It is assumed that the IMU data has AT LEAST the sample frequency of 
% that of the marker data. If the IMU data has a higher sample frequency, 
% some samples will be discarded.
% With option_fs=1, the data is exported at the frame rate of 100Hz by
% default (modify function accordingly for other frame rates)

%% Check input

if ~isfield(Datastr,'IMU') %if no IMU data
    warning(['No field IMU file ' Datastr.Info.Trial '. Skipping.']);
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

% Check if sync indices exist for EMGData
if isfield(Datastr,'EMG')
if ~isfield(Datastr.EMG,'EMGSyncIdx')
    warning('EMG:syncIdx',['No sync indices found in ' Datastr.Info.Trial '. Assuming synchronized EMG data']);
    EMGSyncIdx = [1 size(Datastr.EMG.EMGData,1)];
else
    EMGSyncIdx = Datastr.EMG.EMGSyncIdx;
end
end

% Check if sync indices exist for IMUData
if ~isfield(Datastr.IMU,'IMUSyncIdx')
    warning('IMU:syncIdx',['No sync indices found in ' Datastr.Info.Trial '. Assuming synchronized IMU data']);
    
    IMUSyncIdx = [1 size(Datastr.IMU.IMUData,1)];
else
    IMUSyncIdx = Datastr.IMU.IMUSyncIdx;
end


%%
% Creates motion file of joint angles
% to input in Open Sim from XSens IMU Data


%Depending on option_fs, EMG is exported in different frame rate

disp(option_fs);
switch option_fs
    case 0
nMarkSmpl = markerSyncIdx(end) - markerSyncIdx(1) + 1;
IMUIdx = round(linspace(IMUSyncIdx(1),IMUSyncIdx(end),nMarkSmpl));
IMUdata =  Datastr.IMU.IMUData(IMUIdx,:);
% Frame rates
markerFrameRate = Datastr.Marker.MarkerFrameRate;
IMUtime = (0:size(IMUdata ,1)-1)./markerFrameRate;
allAnglesMatrix=[IMUtime; IMUdata'];
[numRows,numCols] = size(allAnglesMatrix');
ncols = num2str(numCols);
nrows = num2str(numRows);
    case 1
        
% Frame rates
IMUFrameRate = Datastr.IMU.IMUFrameRate;
StandardFrameRate=100; % When using Markerdata, we always export with frameRate=MarkeFrameRate...with the visualeyes this is 100Hz.
[N,D]=rat(StandardFrameRate/IMUFrameRate);

IMUData = resample(Datastr.IMU.IMUData,N,D);
IMUtime = (0:size(IMUData ,1)-1)./StandardFrameRate;%
% 
% figure
% t = (0:size(Datastr.IMU.IMUData ,1)-1)'./IMUFrameRate;%
% plot(t,Datastr.IMU.IMUData(:,10),'+-',IMUtime,IMUData(:,10),'o:')
% legend('original','resampled')


allAnglesMatrix=[IMUtime; IMUData'];
[numRows,numCols] = size(allAnglesMatrix');
ncols = num2str(numCols);
nrows = num2str(numRows);
        
end

% Create base joint angle matrix
%mvnxFile=Datastr.IMU.IMUtrialname;
% splitfilename = strsplit(mvnxFile,'.');
% splitfilename = splitfilename{1};

motFile= [Datastr.IMU.IMUtrialname{1} '.mot'];
mvnxPath = [Datastr.Info.SubjRoot '\' imufolder];
cd(mvnxPath)
fileID = fopen(motFile,'w');
fprintf(fileID,strcat(['first trial\nnRows=',nrows,'\nnColumns=',ncols,'\n\n']));
fprintf(fileID,'# OpenSim Motion File Header:\n');
nrange = num2str([min(min(allAnglesMatrix)) max(max(allAnglesMatrix))]);
otherdata = num2str(1); %always 1 since time vec exists
jointNameStrings=[];
precisionStrings=[];
for i=1:length(Datastr.IMU.IMUDataLabel)
    jointNameStrings=[jointNameStrings '\t' Datastr.IMU.IMUDataLabel{i}];
    precisionStrings=[precisionStrings '%5.4f\t'];
end
jointNameStrings=[jointNameStrings(3:end)];
precisionStrings=[precisionStrings '%5.4f\t\n'];
fprintf(fileID,strcat(['name Joint angles\nnColumns ',ncols,'\nnRows ', nrows,'\ninDegrees=yes','\nendheader\n']));
fprintf(fileID,strcat(['time\t',jointNameStrings,'\n']));
% fprintf(fileID,'%5.3f\t%4.3f\n',Mat);

fprintf(fileID,precisionStrings,allAnglesMatrix);
fprintf(fileID,'\n');
fclose(fileID);
cd(Datastr.Info.SubjRoot)

