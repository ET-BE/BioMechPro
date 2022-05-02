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

function [inData] = getGPi(inData,checkwithFP, folderToSave)
%% Gait Phase Detection based on inertial data
% INPUT)
% inData : - Datastr structure
% 
% OUTPUT)
% If inData is a structure, output is the same structure with added fields
% 
% phasevectorlabel : cell, corresponding with the columns in phasevector
% Contains the name of the phase
% 
% phaseidx : structure with the indices of the start and end of each gait phase
% % 
% The function assumes that the data is smooth and properly filtered !
% 
%
% TODO
%% Store data for output
ContactR = inData.IMU.IMUContactRFoot(:,1) ==1; % Right stance SSR
HSR = find(diff(ContactR)>0)+1;
TOR = find(diff(ContactR)<0);
% Remove ends without start and starts without end
if HSR(end) > TOR(end)
    HSR(end) = [];
end
if TOR(1) < HSR(1)
    TOR(1) = [];
end
peakdist=0.5*inData.IMU.IMUFrameRate;
            valid=(TOR-HSR)>=peakdist;
            HSR=HSR(valid);
            TOR=TOR(valid);


%% Generate output
mvnxPath = [inData.Info.SubjRoot '\' folderToSave];
cd(mvnxPath)
inData.Event.GaitRstanceIMU = [HSR TOR];
csvFile= [inData.IMU.IMUtrialname{1} '.csv'];
csvwrite(csvFile, [inData.IMU.IMUmvnData.time(HSR)', inData.IMU.IMUmvnData.time(TOR)']);
%Plot to check
if checkwithFP
fdata=inData.Force.ForceData;
HSR_FP=inData.Event.GaitRstance(:,1);
TOR_FP=inData.Event.GaitRstance(:,2);
ContactR=inData.IMU.IMUContactRFoot(:,1);
fs_imu=240;
fs_xpc=1000;
figure
plot((1:length(ContactR))/fs_imu*fs_xpc,ContactR*700); hold on
plot(fdata(:,9))
plot(HSR_FP,ones(1,length(HSR_FP))*700,'^')
plot(TOR_FP,ones(1,length(TOR_FP))*700,'v')
plot((HSR)/fs_imu*fs_xpc,ones(1,length(HSR))*700,'o')
plot((TOR)/fs_imu*fs_xpc,ones(1,length(TOR))*700,'*')
legend({'contactR','FzR','HSR FP','TOR FP','HSR IMC','TOR IMC'})
end


end
