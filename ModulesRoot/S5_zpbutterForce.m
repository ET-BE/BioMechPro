function [Datastr] = S5_zpbutterForce(Datastr,fType,fOrd,fCut,keepInput)
% gBMPDynUI fType=1; fOrd=1; fCut=1; keepInput=1;
% 
% Apply butter filter to Force Plates
% 
% INPUT)
% -Datastr: data structure, with at least the field:
% .Force.ForceData
% - fType: string, either 'low','high','bandstop','bandpass'
% 
% - fOrd: integer, filter order
% 
% - fCut: scalar or two element vector, specifying cutoff frequency
% 
% - keepInput: boolean, if true (default) the input will be stored in the
% structure for later use.
%
% OUTPUT)
% -Datastr: data structure, with filtered force data
% 
% NOTES)


%% Check input

if isempty(keepInput)
    keepInput = true;
end

%% Get info
subjtrial = Datastr.Info.Trial;

if ~isfield(Datastr,'Force')
    warning(['Skipping C3D file ' subjtrial ': no force data.']);
    return;
end

fsForce = Datastr.Force.FrameRate;

%% Filter

Datastr.Force.RightForceData(isnan(Datastr.Force.RightForceData(:, 1)), :) = [];  % remove the NaN frames
Datastr.Force.LeftForceData(isnan(Datastr.Force.LeftForceData(:, 1)), :) = [];

[b,a] = butter(fOrd,2.*fCut/fsForce, fType);
Datastr.Force.RightForceData = filtfilt(b,a,Datastr.Force.RightForceData);
Datastr.Force.LeftForceData = filtfilt(b,a,Datastr.Force.LeftForceData);

if keepInput
    Datastr.Info.zpbutterFroce2 = fType;
    Datastr.Info.zpbutterForce3 = fOrd;
    Datastr.Info.zpbutterForce4 = fCut;
end


end
