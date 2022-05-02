function [Datastr]  = I4_FiltIMU(Datastr,fType,fOrd,fCut)
% gBMPDynUI fType=1; fOrd=1; fCut=1;
% 
% INPUT)
% - Datastr: structure, containing at least the fields:
% .Info.SubjRoot
% .IMU.IMUData
% 
% -imuFolder: string, specifying the folder name relative to the root
% folder, that contains the imu measurement for the force data.
% 
% -imuFile: string, specifying a unique identifier of the imu trail, at
% the end of the file.
% 
% - fZeroChan: scalar or vector, containing the channel numbers of the
% force data in the imu measurement, to be subtracted from the force data
% in the data structure.
% 
% - MGRF: matrix, force plate ground reaction force calibration matrix.
% 
% OUTPUT)
% 
% NOTES)
% 


%% Some checks

if ~isfield(Datastr,'IMU')
    warning('No field IMU. Skipping.');
    return;
end

%% Get info

subjrootfolder = Datastr.Info.SubjRoot;
fs_imu = Datastr.IMU.IMUFrameRate;


%% Get filtered imu data

% Get imu
imuData = Datastr.IMU.IMUData;

if ~(isempty(fOrd)||isempty(fCut)||isempty(fType))
    
    % Create filter
    [b,a] = butter(fOrd,2*fCut/fs_imu,fType);

    % Filter imu data
    imuData = filtfilt(b,a,imuData);
    
end

% Store
Datastr.IMU.IMUData = imuData;
 

end