function [Datastr] = S7_zpbutterInsole(Datastr,resampleRate, fType,fOrd,fCut,keepInput)
% gBMPDynUI resampleRate=1; fType=1; fOrd=1; fCut=1; keepInput=1;
% 
% Apply butter filter to Force Plates
% 
% INPUT)
% -Datastr: data structure, with at least the field:
% .Insole.InsoleData
%
% - resampleRate: the resampling rate for the insole data
%
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

if isempty(resampleRate)
    resampleRate = 100;
end

%% Get info
subjtrial = Datastr.Info.Trial;

if ~isfield(Datastr,'Insole')
    warning(['Skipping C3D file ' subjtrial ': no insole data.']);
    return;
end

%% Resampling insole data
Datastr.Insole.ResampleRate = resampleRate;

Datastr.Insole.ResampledLeft = resample(Datastr.Insole.Left(:, 2:end),...
                                Datastr.Insole.Left(:, 1)/1000, resampleRate);
Datastr.Insole.ResampledRight = resample(Datastr.Insole.Right(:, 2:end),...
                                Datastr.Insole.Right(:, 1)/1000, resampleRate);
%% Filter
[b,a] = butter(fOrd,2.*fCut/resampleRate, fType);
Datastr.Insole.FilteredLeft = filtfilt(b,a,Datastr.Insole.ResampledLeft);
Datastr.Insole.FilteredRight = filtfilt(b,a,Datastr.Insole.ResampledRight);

if keepInput
    Datastr.Info.zpbutterInsole2 = fType;
    Datastr.Info.zpbutterInsole3 = fOrd;
    Datastr.Info.zpbutterInsole4 = fCut;
end


end
