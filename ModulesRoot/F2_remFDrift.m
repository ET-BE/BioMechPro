function [Datastr] = F2_remFDrift(Datastr,prctL,prctR)
% gBMPDynUI prctL=1; prctR=1;
% 
% Try to get rid of drift in force channels
% 
% INPUT)
% -Datastr: data structure, with at least the field:
% .Force.ForceData
% containing Nx12 dual plate force data.
% 
% - prctL / prctR: scalar between 1 and 100 (see doc prctile).
% Used to select 'no contact instances' from the vertical forces.
% What is a suitable prctile depends your data distribution, which in turn
% depends on the gait frequency. Suggestion: something between 5 and 25.
% 
% OUTPUT)
% -Datastr: data structure, with detrended force data
% 
% NOTES)
% Input is assumed an Nx12 matrix with dual plate force data, and the
% vertical GRF data being in channels 3 (left) and 9 (right).

%% Remove drift

% Remove drift and store
if isfield(Datastr,'Force')
    Datastr.Force.ForceData = remFDrift(Datastr.Force.ForceData,prctL,prctR);
else
    warning('F2_remFDrift:nofield','No field Force. Skipping.');
end

end