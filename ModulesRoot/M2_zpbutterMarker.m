function [Datastr] = M2_zpbutterMarker(Datastr,fType,fOrd,fCut,keepInput)
% gBMPDynUI fType=1; fOrd=1; fCut=1; keepInput=1;
% 
% Filter marker data with a zero phase butterworth filter
% 
% INPUT)
% - Datastr, structure with the field:
% .Marker.MarkerData
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
% - Datastr, structure with added fields (if keepInput = true):
% .Info.zpbutterMarker2 (fType)
% .Info.zpbutterMarker3 (fOrd)
% .Info.zpbutterMarker4 (fCut)
% 
% NOTES)

%% Check input

if isempty(keepInput)
    keepInput = true;
end

%% Get info
subjtrial = Datastr.Info.Trial;

if ~isfield(Datastr,'Marker')
    warning(['Skipping C3D file ' subjtrial ': no marker data.']);
    return;
end

fsMark = Datastr.Marker.MarkerFrameRate;

%% Filter

[b,a] = butter(fOrd,2.*fCut/fsMark,fType);
Datastr.Marker.MarkerData = filtfilt(b,a,Datastr.Marker.MarkerData);

if keepInput
    Datastr.Info.zpbutterMarker2 = fType;
    Datastr.Info.zpbutterMarker3 = fOrd;
    Datastr.Info.zpbutterMarker4 = fCut;
end

end