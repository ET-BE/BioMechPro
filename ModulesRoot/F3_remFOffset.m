function [Datastr] = F3_remFOffset(Datastr,offset)
% gBMPDynUI offset=1;
% 
% Specify offset of Force Plates
% 
% INPUT)
% -Datastr: data structure, with at least the field:
% containing Nx12 dual plate force data. e.g.
% offset=[0 -0.1500 0 0 -0.1500 0]; In OS coordinates z x y. This example
% would subtract 15 cm in the x direction (forward walking) in both plates
% 
% OUTPUT)
% -Datastr: data structure, with detrended force data
% 
% NOTES)


%% Remove offset

% Remove drift and store
if isfield(Datastr,'Force')
    Datastr.Force.ForcePlateOffset=offset;
else
    warning(['No field Force in trial ' Datastr.Info.Trial '. Skipping.']);
    return
end
