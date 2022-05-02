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

function [Datastr] = S2_importIMU(Datastr,imufolder,frommvnx, removemvnx,removemat,trunkmotion, translation)
% gBMPDynUI imufolder=1; frommvnx=1; removemvnx=1; removemat=1; trunkmotion=1; translation=1;
% Create a .mot file with join angles for OpenSim input from Xsens data (.mvnx)
% INPUT)
% imufolder=['IMU']; specifying the folder containing the measurement
% trials with the IMU data, relative to the subject root folder.
% frommvnx=1; %Start from mvnx file (if 0 start from mat file)
% removemvnx=0; %Free some space as .mvnx files are big
% removemat=0; %Free some space as .mat files are big
% translation=0; %=0 no translation of pelvis, =1 with no pelvis_tx, pelvis_ty or pelvis_tz
% trunkmotion=0; %=0 nothing, =1 with jL5S1, =2 with quaternions (set values for...
% pelvis_tilt, pelvis_list, pelvis_rotation and lumbar_extension,lumbar_bending,lumbar_rotation)

%
% -removemvnx: boolean, if true it removes .mvnx files in folder to free some space
% -trunkmotion: boolean, get values for pelvis and lumbar or set them to 0
% if false
% OUTPUT)
% - Datastr, structure, with added fields:
% .IMU.IMUData
% .IMU.IMUDataLabel
% .IMU.IMUFrameRate
%
% NOTES)
%

%% Check input
%Datastr.Info.SubjRoot %Datastr.Info.Trial
mvnxPath = [Datastr.Info.SubjRoot '\' imufolder];
cd(mvnxPath)

if frommvnx
    mvnxfiles = dir('*.mvnx');
    if isempty(mvnxfiles)
        warning(['Skipping IMU processing, *.mvnx files not found in IMU folder:' mvnxPath]);
        return
    end
    for ifile = 1:size(mvnxfiles,1)
        mvnxfilenames{ifile} = mvnxfiles(ifile).name;
        % Reverse string and remove '.mvnx', for comparison with strncmpi
        mvnxfilenamescmp{ifile} = mvnxfiles(ifile).name(end-5:-1:1);        
    end
    clear mvnxfiles;
    match = strncmpi( Datastr.Info.Trial(end:-1:1),mvnxfilenamescmp,length(Datastr.Info.Trial) );
    if sum(match) == 1
        mvnxFile = mvnxfilenames{match};
        [mvnData] = loadMVNX([mvnxFile]); %Modified function loadMVNXdata.m from Xsens MVN Studio Developer Toolkit 1.0.2. Contains load_mvnx.m function from same kit.
        split=strsplit(mvnxFile,'.');
        trialname=split(1);
    else
        warning(['Skipping file ' Datastr.Info.Trial ': not found or m ultiple found.']);
        return;
    end
    if ~removemat
        matFile = strrep(mvnxFile,'.mvnx','.mat'); %Keep same name
        save(matFile,'mvnData');
    end
    
else
    matfiles = dir('*.mat');
    if isempty(matfiles)
        warning(['Skipping IMU processing, *.mat files not found in IMU folder:' mvnxPath]);
        return
    end
    for ifile = 1:length(matfiles)
        matfilenames{ifile} = matfiles(ifile).name;
        % Reverse string and remove '.mat', for comparison with strncmpi
        matfilenamescmp{ifile} = matfiles(ifile).name(end-4:-1:1);
    end
    clear matfiles;
    match = strncmpi( Datastr.Info.Trial(end:-1:1),matfilenamescmp,length(Datastr.Info.Trial) );
    if sum(match) == 1
        matFile = matfilenames{match};
        %matFile = strrep(mvnxFile,'.mvnx','.mat'); %Keep same name
        load(matFile) %loads MVNData
        split=strsplit(matFile,'.');
        trialname=split(1);
    else
        warning(['Skipping file ' Datastr.Info.Trial ': not found or multiple found.']);
        return;
    end
    
end




% Process current trial


% Based on read_jointAngles_lowerbody_write_motion_file.m script (found online)
%% Create table with jointvalues
jointLabelsMVN = fieldnames(mvnData.joints.jointAngle);
[nJoints,~] = size(jointLabelsMVN);
jointNames = cell(3*nJoints+1,1);

for hh = 1:nJoints
    jointName = jointLabelsMVN{hh};
    %repeat jointName 3 times
    jointNames{hh*3-2,1} = jointName;
    jointNames{hh*3-1,1} = jointName;
    jointNames{hh*3,1} = jointName;
    %add time if you're at the last joint
    if hh == nJoints
        jointNames{hh*3+1,1} = cellstr('time');
    end
end

ColIndx = repmat([1,2,3]',nJoints,1); %These correspond to X,Y,Z angles (ColIndx = 1,2,3)
ColIndx = [ColIndx; 0]; %add a zero for time

%Bad way of gettting numel of a joint angle
[nSamples,~] = size(mvnData.joints.jointAngle.(jointLabelsMVN{1}));

Angles = zeros(length(jointNames),nSamples);  %Angles in degrees
jointTable = table(jointNames,ColIndx,Angles);

axis=0;
jj=1;
for ii = 1:length(jointTable.jointNames)-1
    axis=axis+1;
    switch axis
        %     for sampleIndex = 4:nSamples
        %         temp(:,sampleIndex)= tree.subject.frames.frame(sampleIndex).jointAngle(ii); %joint angles
        %     end
        case 1
            jointTable.Angles(ii,:) = mvnData.joints.jointAngle.(jointLabelsMVN{jj})(:,1);%X
        case 2
            jointTable.Angles(ii,:) = mvnData.joints.jointAngle.(jointLabelsMVN{jj})(:,2);%Y
        case 3
            jointTable.Angles(ii,:) = mvnData.joints.jointAngle.(jointLabelsMVN{jj})(:,3);%Z
            axis=0;
            jj=jj+1;
    end
end
%Create time vector and add to table
% Time = zeros(nSamples,1);
% for jj=[1:nSamples]
%     Time(jj)=tree.subject.frames.frame(jj).time;
% end
% Time=Time./1000;
% fs=240; %sampling frequency from MVN link
Time=mvnData.time;
jointTable.Angles(length(jointNames),:) = Time;

%This gives you a joint table with each joint segment (as labeled in MVN
%file and then gives you a row vector of each X,Y,Z angles for each joint over time)
% Now you can pull out vectors that you want


%% Pull out the angles we want for this model
% E.g.
% pelvis_tilt	= Z (+/-)
% pelvis_tx	= x (position 1 x)
% pelvis_ty	= y  (position 2 y)
% hip_flexion_r	= Z (+/-)
% knee_angle_r = only Z
% ankle_angle_r = Z

%% first get the orientation and vertical translation of the pelvis
[z,y,x]=quat2angle(mvnData.segments.orientation.Pelvis(:,:));

% if the orentations are close to pi/-pi, then compensate them to zero
for io = 1:length(x)
%     if x(io) > pi/2
%         x(io) = x(io) - pi;
%     elseif x(io) < -pi/2
%         x(io) = x(io) + pi;
%     end
%     
%     if y(io) > pi/2
%         y(io) = y(io) - pi;
%     elseif y(io) < -pi/2
%         y(io) = y(io) + pi;
%     end
    
    if z(io) < 0
        z(io) = z(io) + 2*pi;
    end
end

% assign the extracted & compenstated pelvis orientations
    
pelvis_tilt = -y'*180/pi; %in Xsens, the positive rotation is around the oppsite direction of Y
pelvis_list = x'*180/pi;
pelvis_rotation = z'*180/pi;
pelvis_tx = zeros(1,nSamples); % -(mvnData.segments.position.Pelvis(:,1) - mvnData.segments.position.Pelvis(1,1))';  % Xense avatar is walking the opposite direction of the X
pelvis_ty = mvnData.segments.position.Pelvis(:,3)';  
pelvis_tz = zeros(1,nSamples); % (mvnData.segments.position.Pelvis(:,2) - mvnData.segments.position.Pelvis(1,2))'; % Xense avatar is drifting the left side

%% then get the trunk joint angles

negZ = -1;

%Method 1) adding L5S1, L4L3, L1T12, T9T8 joint angle together
rows1 = strcmp(jointTable.jointNames,'jL5S1') & jointTable.ColIndx==3;
rows2 = strcmp(jointTable.jointNames,'jL4L3') & jointTable.ColIndx==3;
rows3 = strcmp(jointTable.jointNames,'jL1T12') & jointTable.ColIndx==3;
rows4 = strcmp(jointTable.jointNames,'jT9T8') & jointTable.ColIndx==3;
vars = {'Angles'};
% in OpenSim the extension is positive, then a negative sign is added here
lumbar_extension = -(jointTable{rows1,vars} + jointTable{rows2,vars} +...
                  jointTable{rows3,vars} + jointTable{rows4,vars});


rows1 = strcmp(jointTable.jointNames,'jL5S1') & jointTable.ColIndx==1;
rows2 = strcmp(jointTable.jointNames,'jL4L3') & jointTable.ColIndx==1;
rows3 = strcmp(jointTable.jointNames,'jL1T12') & jointTable.ColIndx==1;
rows4 = strcmp(jointTable.jointNames,'jT9T8') & jointTable.ColIndx==1;
vars = {'Angles'};   

lumbar_bending = jointTable{rows1,vars} + jointTable{rows2,vars} +...
                  jointTable{rows3,vars} + jointTable{rows4,vars};
rows1 = strcmp(jointTable.jointNames,'jL5S1') & jointTable.ColIndx==2;
rows2 = strcmp(jointTable.jointNames,'jL4L3') & jointTable.ColIndx==2;
rows3 = strcmp(jointTable.jointNames,'jL1T12') & jointTable.ColIndx==2;
rows4 = strcmp(jointTable.jointNames,'jT9T8') & jointTable.ColIndx==2;
vars = {'Angles'};
lumbar_rotation = jointTable{rows1,vars} + jointTable{rows2,vars} +...
                  jointTable{rows3,vars} + jointTable{rows4,vars};



%Rest of vars
rows = strcmp(jointTable.jointNames,'jRightHip') & jointTable.ColIndx==3;
vars = {'Angles'};
hip_flexion_r = jointTable{rows,vars};

rows = strcmp(jointTable.jointNames,'jRightHip') & jointTable.ColIndx==1;
vars = {'Angles'};
% in xsens, abduction is positive
hip_adduction_r = -jointTable{rows,vars};

rows = strcmp(jointTable.jointNames,'jRightHip') & jointTable.ColIndx==2;
vars = {'Angles'};
hip_rotation_r = jointTable{rows,vars};

rows = strcmp(jointTable.jointNames,'jRightKnee') & jointTable.ColIndx==3;
vars = {'Angles'};
% in opensim, keen flexion is negative
knee_angle_r = negZ*jointTable{rows,vars};

%We dont use this in our model in OpenSim
rows = strcmp(jointTable.jointNames,'jRightKnee') & jointTable.ColIndx==1;
vars = {'Angles'};
% in xsens abduction is positive.
knee_adduction_r = negZ*jointTable{rows,vars};

%We dont use this in our model in OpenSim
rows = strcmp(jointTable.jointNames,'jRightKnee') & jointTable.ColIndx==2;
vars = {'Angles'};
% internal rotation is positive
knee_rotation_r = jointTable{rows,vars};

rows = strcmp(jointTable.jointNames,'jRightAnkle') & jointTable.ColIndx==3;
vars = {'Angles'};
% Both in OpenSim and Xsens, Dosiflexion is positive
ankle_angle_r = jointTable{rows,vars};

%We dont use this in our model in OpenSim
rows = strcmp(jointTable.jointNames,'jRightAnkle') & jointTable.ColIndx==1;
vars = {'Angles'};
% in opensim, adduction is positive
ankle_adduction_r = negZ*jointTable{rows,vars};

rows = strcmp(jointTable.jointNames,'jRightAnkle') & jointTable.ColIndx==2;
vars = {'Angles'};
subtalar_angle_r = jointTable{rows,vars};
%Checking ik.mot from camera this angle is set to 0. Keep same for IMU!
% rows = strcmp(jointTable.jointNames,'jRightBallFoot') & jointTable.ColIndx==3;
% vars = {'Angles'};
% mtp_angle_r = jointTable{rows,vars};
mtp_angle_r =zeros(1,nSamples);

rows = strcmp(jointTable.jointNames,'jLeftHip') & jointTable.ColIndx==3;
vars = {'Angles'};
hip_flexion_l = jointTable{rows,vars};

rows = strcmp(jointTable.jointNames,'jLeftHip') & jointTable.ColIndx==1;
vars = {'Angles'};
hip_adduction_l = negZ*jointTable{rows,vars};

rows = strcmp(jointTable.jointNames,'jLeftHip') & jointTable.ColIndx==2;
vars = {'Angles'};
hip_rotation_l = jointTable{rows,vars};

rows = strcmp(jointTable.jointNames,'jLeftKnee') & jointTable.ColIndx==3;
vars = {'Angles'};
knee_angle_l = negZ*jointTable{rows,vars};
%We dont use this in our model in OpenSim
rows = strcmp(jointTable.jointNames,'jLeftKnee') & jointTable.ColIndx==1;
vars = {'Angles'};
knee_adduction_l = negZ*jointTable{rows,vars};
%We dont use this in our model in OpenSim
rows = strcmp(jointTable.jointNames,'jLeftKnee') & jointTable.ColIndx==2;
vars = {'Angles'};
knee_rotation_l = jointTable{rows,vars};

rows = strcmp(jointTable.jointNames,'jLeftAnkle') & jointTable.ColIndx==3;
vars = {'Angles'};
ankle_angle_l = jointTable{rows,vars};
%We dont use this in our model in OpenSim
rows = strcmp(jointTable.jointNames,'jLeftAnkle') & jointTable.ColIndx==1;
vars = {'Angles'};
ankle_adduction_l = negZ*jointTable{rows,vars};

rows = strcmp(jointTable.jointNames,'jLeftAnkle') & jointTable.ColIndx==2;
vars = {'Angles'};
subtalar_angle_l = jointTable{rows,vars};

%Checking ik.mot from camera this angle is set to 0. Keep same for IMU!
% rows = strcmp(jointTable.jointNames,'jLeftBallFoot') & jointTable.ColIndx==3;
% vars = {'Angles'};
% mtp_angle_l = jointTable{rows,vars};
mtp_angle_l =zeros(1,nSamples);

%% Extract the CoM
CoMx = -(mvnData.centerOfMass(:, 1) - mvnData.centerOfMass(1, 1));
CoMy = mvnData.centerOfMass(:, 3);
CoMz = (mvnData.centerOfMass(:, 2) - mvnData.centerOfMass(1, 2));
%%
%Now output all to matrix

if isfield(mvnData,'contact')
    Datastr.IMU.IMUContactRFoot = mvnData.contact;
end

Datastr.IMU.IMUData = [pelvis_tilt; pelvis_list; pelvis_rotation; pelvis_tx;  pelvis_ty; pelvis_tz; hip_flexion_r; hip_adduction_r; hip_rotation_r; knee_angle_r; knee_adduction_r; knee_rotation_r; ankle_angle_r; ankle_adduction_r; subtalar_angle_r; mtp_angle_r; hip_flexion_l; hip_adduction_l; hip_rotation_l; knee_angle_l; knee_adduction_l; knee_rotation_l; ankle_angle_l; ankle_adduction_l; subtalar_angle_l; mtp_angle_l; lumbar_extension; lumbar_bending; lumbar_rotation]';
Datastr.IMU.IMUDataLabel = {'pelvis_tilt' 'pelvis_list' 'pelvis_rotation' 'pelvis_tx' 'pelvis_ty' 'pelvis_tz' 'hip_flexion_r' 'hip_adduction_r' 'hip_rotation_r' 'knee_angle_r' 'knee_adduction_r' 'knee_rotation_r' 'ankle_angle_r' 'ankle_adduction_r' 'subtalar_angle_r' 'mtp_angle_r' 'hip_flexion_l' 'hip_adduction_l' 'hip_rotation_l' 'knee_angle_l' 'knee_adduction_l' 'knee_rotation_l' 'ankle_angle_l' 'ankle_adduction_l' 'subtalar_angle_l' 'mtp_angle_l' 'lumbar_extension' 'lumbar_bending' 'lumbar_rotation'};
Datastr.IMU.CoM = [CoMx, CoMy, CoMz];
Datastr.IMU.Sensor = mvnData.sensors;
Datastr.IMU.IMUFrameRate = round( 1 ./ mean(diff(mvnData.time(1:end-1))) );%% should be  sampling frequency from MVN link (240Hz)
%only if modified loadMVNXdata.m: Datastr.IMU.IMUFrameRate = mvnData.frameRate

% Datastr.IMU.IMUmvnData = mvnData; %for debugging, raw matrix
Datastr.IMU.IMUtrialname = trialname;
Datastr.IMU = orderfields(Datastr.IMU);

if removemvnx %Free some space by removing .mvnx files
    delete(mvnxFile)
end
if removemat
    delete(matFile)
end
cd(Datastr.Info.SubjRoot)
end