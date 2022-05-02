function [Datastr] = S4_zpbutterMarker(Datastr,fType,fOrd,fCut,keepInput)
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

fsMark = Datastr.Marker.FrameRate;

%% Filter

[b,a] = butter(fOrd,2.*fCut/fsMark,fType);
for mk = 1:length(Datastr.Marker.MarkerData(:, 1, 1))
    
    MrData = squeeze(Datastr.Marker.MarkerData(mk, :, :));  % extract each marker data
    
    id = find(isnan(MrData(1, :))==0);  % find ~NaN data frames
    id_nan = find(isnan(MrData(1, :))==1);  % find NaN data frames
    brk_id = find(diff(id) > 1);  % find the breaking points of ~NaN data periods
    
    % filt each data period
    if ~isempty(brk_id)
        for bp = 1:length(brk_id)+1
            if bp == 1
                if id(brk_id(bp)) - id(1) < 6
                    continue
                else
                Datastr.Marker.MarkerData(mk, :, id(1):id(brk_id(bp))) =...
                    filtfilt(b,a,MrData(:, id(1):id(brk_id(bp)))')';
                end
            end
            if bp == length(brk_id)
                if id(end) - id(brk_id(bp)+1) < 6
                    continue
                else
                Datastr.Marker.MarkerData(mk, :, id(brk_id(bp)+1):id(end)) =...
                    filtfilt(b,a,MrData(:, id(brk_id(bp)+1):id(end))')';
                end
            end
            if bp > 1 && bp < length(brk_id)
                if id(brk_id(bp+1)) - id(brk_id(bp)+1) < 6
                    continue
                else
                Datastr.Marker.MarkerData(mk, :, id(brk_id(bp)+1):id(brk_id(bp+1))) =...
                    filtfilt(b,a,MrData(:, id(brk_id(bp)+1):id(brk_id(bp+1)))')';
                end
            end
        end
    end

    Datastr.Marker.MarkerData(mk, :, id_nan) = NaN;  % assign NaN to the origonal NaN frames

end

if keepInput
    Datastr.Info.zpbutterMarker2 = fType;
    Datastr.Info.zpbutterMarker3 = fOrd;
    Datastr.Info.zpbutterMarker4 = fCut;
end

end