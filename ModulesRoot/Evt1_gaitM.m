function [Datastr] = Evt1_gaitM(Datastr,feetCell)
% gBMPDynUI feetCell=1;
% 
% See getGPm for help file
% There are some assumptions on the Marker data dimensions
% 
% Example value for feetCell:
% {'CAL','CM1L','CAR','CM1R'}

Datastr = getGPm(Datastr,feetCell);

end