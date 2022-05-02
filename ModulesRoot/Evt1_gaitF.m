function [Datastr] = Evt1_gaitF(Datastr,fthresh)
% gBMPDynUI fThresh=1;
% 
% See getGPf for help file
% Works for dual force plate only

% Get gait information based on forces
if isfield(Datastr,'Force')
    Datastr = getGPf(Datastr,fthresh);
else
    warning(['No field Force in trial ' Datastr.Info.Trial '. Skipping']);
end

end