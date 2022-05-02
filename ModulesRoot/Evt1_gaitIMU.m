function [Datastr] = Evt1_gaitIMU(Datastr,checkwithFP, folderToSave)
% gBMPDynUI checkwithFP=1; folderToSave=1;
%% Gait phase detection based on IMU right foot contact data
% Get gait information based on IMU
if isfield(Datastr,'IMU')
    Datastr = getGPi(Datastr,checkwithFP, folderToSave);
else
    warning(['No field IMU in trial ' Datastr.Info.Trial '. Skipping']);
end

end
