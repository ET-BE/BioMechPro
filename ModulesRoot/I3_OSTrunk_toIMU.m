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

function [Datastr] = I3_OSTrunk_toIMU(Datastr,osFolder)
% gBMPDynUI osFolder=1;
% oFolder=['OS];
% Include  trunk values in IMU data from IK of camera in OpenSim
% INPUT)
% - osFolder: string, specifying the folder containing the measurement
% trials with the IMU data, relative to the subject root folder.
% Example:
% 'OS'
%
% OUTPUT)
% - Datastr, structure, with added fields:
% .IMU.IMUData
%
%
% This function subtitutes IMU data values by the data values obtained through marker and OpenSim for the folliwng fields:
% 'pelvis_tx','pelvis_ty','pelvis_tz','pelvis_list','pelvis_rotation','pelvis_tilt','lumbar_bending','lumbar_rotation','lumbar_extension' 
% This allows IMU data to have same origin and trunk motion as the one
% derived from the optical system. Important if we want to include force
% information toinertial data with same coordinate frame.

% Needs further refining


%% Check
if ~isfield(Datastr,'IMU') % if no IMU data
    warning(['No IMU file ' Datastr.Info.Trial '. Skipping.']);
    return
end

if ~isfield(Datastr,'Marker') % if no VZ data
    warning(['No VZ data in file ' Datastr.Info.Trial '. Skipping.'])
    return
end

% Check if sync indices exist for MarkerData and IMUData
if ~isfield(Datastr.Marker,'MarkerSyncIdx') || ~isfield(Datastr.IMU,'IMUSyncIdx')
    warning('IMU:syncIdx',['No sync indices found in ' Datastr.Info.Trial '. Synchronize first marker-IMU data'])
    markerSyncIdx = [1 size(Datastr.Marker.MarkerData,1)];
    IMUSyncIdx = [1 size(Datastr.IMU.IMUData,1)];
else
    markerSyncIdx = Datastr.Marker.MarkerSyncIdx;
    IMUSyncIdx = Datastr.IMU.IMUSyncIdx;
end

%% Get info and paths
rootfolder = Datastr.Info.SubjRoot;
bckslsh = strfind(rootfolder,'\');
if isempty(bckslsh)
    bckslsh = 0;
end

trial = Datastr.Info.Trial;
savename = [rootfolder(bckslsh(end)+1:end) trial];
% Get path to created IK.mot file
subjosikmot = [rootfolder '\' osFolder '\' savename 'IK.mot'];

if ~isfield(Datastr.Marker,'JointAngData') %Check IK was performed with camera system in OpenSim. Assuming if field JointAngData exist,the correspoding .mot file exists as well
    warning(['No .mot (from OpenSim) data in file ' Datastr.Info.Trial '. Skipping.'])
    return
end


%% Load data

% Read file header and data header
fid = fopen(subjosikmot,'r');
if fid == -1
    error('getOSIKinM:fileopen',['Failed to open ' subjosikmot]);
end


tline = '';
ntimeout = 0;
while isempty(strfind(tline,'time')) && ntimeout < 20
    
    tline = fgetl(fid);
    
    if ~isempty(strfind(tline,'nRows='))
        nRows = str2double( tline(strfind(tline,'nRows=') + 6 : end) );
    elseif ~isempty(strfind(tline,'nColumns='))
        nCols = str2double( tline(strfind(tline,'nColumns=') + 9 : end) );
    elseif ~isempty(strfind(tline,'inDegrees='))
        isDeg = tline(strfind(tline,'inDegrees=') + 10 : end);
    end
    
    ntimeout = ntimeout + 1;
end
fclose(fid);

if ntimeout == 20
    error('getOSIKinM:DataHeader','Unable to find data column header in file.');
else
    angData = dlmread(subjosikmot,'\t',[ntimeout,0,ntimeout+nRows-1,nCols-1]);
    % NOTE: assumed here the data is directly below the column header
end

%% Sort data based on header
%
% NOTE: It is assumed the column header is delimited with tabs
% NOTE2: Even though the torques are not in global coordinates, they are
% sorted in a "global as possible" way. For example, the knee
% flexion-extension moment is assumed to be in the sagittal plane, about
% the transversal axis, regardless of e.g. the upper leg rotation.
% As a reminder, in OpenSim, x is the walking direction, y is upward and z
% to the right

headercell = strsplit(tline,'\t');

% Introduce angData at the right times in IMUdata
%(i.e. take into account that IMU and VZ data have diff sampling rate and therefore values from the .mot must be inserted at specific times in IMU data)
% nMarkSmpl = length(angData);
nMarkSmpl = markerSyncIdx(end) - markerSyncIdx(1) + 1;
IMUIdx = round(linspace(IMUSyncIdx(1),IMUSyncIdx(end),nMarkSmpl));

% Create indices and labelcells for data storage based on header names
for ihd = 1:length(headercell)
    switch headercell{ihd}
        case {'pelvis_tx','pelvis_ty','pelvis_tz','pelvis_list','pelvis_rotation','pelvis_tilt','lumbar_bending','lumbar_rotation','lumbar_extension'}         
            iDx= strcmp(Datastr.IMU.IMUDataLabel,headercell{ihd});
                Datastr.IMU.IMUData(IMUIdx(1:length(angData)),iDx)=angData(:,ihd);
        otherwise % time / other joints
            % Do nothing
    end
end
