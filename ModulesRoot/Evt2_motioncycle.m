function [Datastr] = Evt2_motioncycle(Datastr,Datafrom)
% gBMPDynUI Datafrom=1;
%
% Detect start of each repetition of the following motor tasks: gait,
% squads and calf raise. 

% Datafrom= The detection is based on the the joint angles
% derived through the marker data or the IMU. 
% Example:
% Datafrom= 'Marker'


% Get motion cycle information based on Marker
if strcmp(Datafrom,'Marker')
    if isfield(Datastr,'Marker')
        Datastr = getmotioncycle(Datastr,Datafrom);
        Datastr.Event.MotionCyclefrom = Datafrom;
    else
        warning(['No field Marker in trial ' Datastr.Info.Trial '. Skipping']);
    end
    
    
elseif strcmp(Datafrom,'IMU')
    if isfield(Datastr,'IMU')
        Datastr = getmotioncycle(Datastr,Datafrom);
        Datastr.Event.MotionCyclefrom = Datafrom;
    else
        warning(['No field IMU in trial ' Datastr.Info.Trial '. Skipping']);
    end
end
end