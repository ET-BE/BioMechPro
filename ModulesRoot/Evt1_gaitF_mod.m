function [Datastr] = Evt1_gaitF_mod(Datastr,fthresh)
% gBMPDynUI fThresh=1;
% 

% Modfied version of Evt1_gaitF: prevents crashing with no detection of
% gait events. For example trials of squats or calf raise

% See getGPf for help file
% Works for dual force plate only

% Get gait information based on forces
if isfield(Datastr,'Force')
    Datastr = getGPf_mod(Datastr,fthresh);
else
    warning(['No field Force in trial ' Datastr.Info.Trial '. Skipping']);
end

end