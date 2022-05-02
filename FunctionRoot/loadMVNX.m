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

function [mvnData] = loadMVNX(mvnxPath)
% addpath('C:\Program Files\Xsens\Xsens MVN Developer Toolkit 1.0.7\Matlab\')
%mvnxPath = 'D:\KNEEMO Dataset\VUMC Dataset\2017_07\MVN\HoloLens-014.mvnx';
mvnx = load_mvnx(mvnxPath);
%% Get segment labels
segLabels = {mvnx.subject.segments.segment.label}';
if isfield(mvnx.subject, 'sensors' )
    senLabels = {mvnx.subject.sensors.sensor.label}';
else
    senLabels = {};
end
if isfield(mvnx.subject, 'joints' )
    jntLabels = {mvnx.subject.joints.joint.label}';
else
    jntLabels = {};
end

% Get number of segments
segNum = size(segLabels,1);
senNum = size(senLabels,1);
jntNum = size(jntLabels,1);
% Get time including n-pose, identity pose and t-pose data
timeNum = size(mvnx.subject.frames.frame,2);
% Loop for time and for segments
for iTime = 4:timeNum %skips calibration data
    iDx=iTime-3; %Start with index 1 (when iTime starts with 4 (iTime-3)=1)
    % Time Data %
    %Added_Laura
    mvnData.time(iDx)=mvnx.subject.frames.frame(iTime).time; 
    %
    % Segment Data %
    for iSeg = 1:segNum
        % Angular Variables
        iVar = (iSeg-1)*4+1; % Nx4 for quaternion orientation
        if isfield(mvnx.subject.frames.frame, 'orientation')
            mvnData.segments.orientation.(segLabels{iSeg})(iDx,1:4) = mvnx.subject.frames.frame(iTime).orientation(iVar:iVar+3);
        end
        iVar = (iSeg-1)*3+1; % Nx3
        if isfield(mvnx.subject.frames.frame, 'angularVelocity')
            mvnData.segments.angularVelocity.(segLabels{iSeg})(iDx,1:3) = mvnx.subject.frames.frame(iTime).angularVelocity(iVar:iVar+2);
        end
        if isfield(mvnx.subject.frames.frame, 'angularAcceleration')
            mvnData.segments.angularAcceleration.(segLabels{iSeg})(iDx,1:3) = mvnx.subject.frames.frame(iTime).angularAcceleration(iVar:iVar+2);
        end
        % Linear Variables
        if isfield(mvnx.subject.frames.frame, 'position')
            mvnData.segments.position.(segLabels{iSeg})(iDx,1:3) = mvnx.subject.frames.frame(iTime).position(iVar:iVar+2);
        end
        if isfield(mvnx.subject.frames.frame, 'velocity')
            mvnData.segments.velocity.(segLabels{iSeg})(iDx,1:3) = mvnx.subject.frames.frame(iTime).velocity(iVar:iVar+2);
        end
        if isfield(mvnx.subject.frames.frame, 'acceleration')
            mvnData.segments.acceleration.(segLabels{iSeg})(iDx,1:3) = mvnx.subject.frames.frame(iTime).acceleration(iVar:iVar+2);
        end
%         if isfield(mvnx.subject.frames.frame, 'acceleration')
%             mvnData.segments.acceleration.(segLabels{iSeg})(iDx,1:3) = mvnx.subject.frames.frame(iTime).acceleration(iVar:iVar+2);
%         end
    end
    % Sensor Data %
    for iSen = 1:senNum
        iVar = (iSen-1)*4+1; % Nx4
        if isfield(mvnx.subject.frames.frame, 'sensorOrientation')
            mvnData.sensors.sensorOrientation.(senLabels{iSen})(iDx,1:4) = mvnx.subject.frames.frame(iTime).sensorOrientation(iVar:iVar+3);
        end
        iVar = (iSen-1)*3+1; % Nx3
        if isfield(mvnx.subject.frames.frame, 'sensorAcceleration')
            mvnData.sensors.sensorAcceleration.(senLabels{iSen})(iDx,1:3) = mvnx.subject.frames.frame(iTime).sensorAcceleration(iVar:iVar+2);
        end
        if isfield(mvnx.subject.frames.frame, 'sensorAngularVelocity')
            mvnData.sensors.sensorAngularVelocity.(senLabels{iSen})(iDx,1:3) = mvnx.subject.frames.frame(iTime).sensorAngularVelocity(iVar:iVar+2);
        end
        if isfield(mvnx.subject.frames.frame, 'sensorMagneticField')
            mvnData.sensors.sensorMagneticField.(senLabels{iSen})(iDx,1:3) = mvnx.subject.frames.frame(iTime).sensorMagneticField(iVar:iVar+2);
        end
    end
    % Joint Data %
    for iJnt = 1:jntNum
        iVar = (iJnt-1)*3+1; % Nx3
        if isfield(mvnx.subject.frames.frame, 'jointAngle')
            mvnData.joints.jointAngle.(jntLabels{iJnt})(iDx,1:3) = mvnx.subject.frames.frame(iTime).jointAngle(iVar:iVar+2);
        end
        if isfield(mvnx.subject.frames.frame, 'jointAngleXZY')
            mvnData.joints.jointAngleXZY.(jntLabels{iJnt})(iDx,1:3) = mvnx.subject.frames.frame(iTime).jointAngleXZY(iVar:iVar+2);
        end
    end
    if isfield(mvnx.subject.frames.frame, 'centerOfMass')
        mvnData.centerOfMass(iDx,1:3) = mvnx.subject.frames.frame(iTime).centerOfMass(1:3);
    end
    %Added_Laura: contact
    if isfield(mvnx.subject.frames.frame,'contact')
        if ~iscell(mvnx.subject.frames.frame(iTime).contact)
            mvnData.contact(iDx,1)=NaN;
        else
            if strcmp('RightFoot',mvnx.subject.frames.frame(iTime).contact{1,1})
                mvnData.contact(iDx,1)=1;
            elseif strcmp('LeftFoot',mvnx.subject.frames.frame(iTime).contact{1,1})
                mvnData.contact(iDx,1)=0;
            end
        end
    end
    %
end
%Added_Laura
mvnData.time=mvnData.time./1000; %to seconds
%
mvnData.points = mvnx.subject.segments.segment;
end %