function [Datastr] = A1_sortAnalogC3D(Datastr,fChan,fLabel,eChan,eLabel,oChan,oLabel)
% gBMPDynUI forceChan=1; forceLbl=1; emgChan=1; emgLbl=1; otherChan=1; otherLbl=1;
%
% Sort analog data channels into corresponding fields of the data structure
% 
% INPUT)
% - Datastr, structure, with at least the field:
% .Analog.AnalogData
% 
% - fChan, eChan, oChan, scalars or vectors, containing the indices
% (channel numbers) in which force, EMG, and any other data are stored
% respectively. Each field can be left empty if none are available.
% 
% - fLabel, eLabel, oLabel, cells, containing labels with names of the
% channels. (TODO: these can also be extracted from the C3D if channel
% names were given in the data acquisition software)
% 
% OUTPUT)
% - Datastr, structure, with one or more of the following added fields:
% .Force.ForceData
% .Force.ForceDataLabel
% .Force.ForceFrameRate
% .EMG.EMGData
% .EMG.EMGDataLabel
% .EMG.EMGFrameRate
% .Other.OtherData
% .Other.OtherDataLabel
% .Other.OtherFrameRate
% 
% NOTES)
% 

%% Sort the data

% Force data
if ~isempty(fChan)
    
    Datastr.Force.ForceData = Datastr.Analog.AnalogData(:,fChan);
    Datastr.Force.ForceFrameRate = Datastr.Analog.AnalogFrameRate;
    
    %Datastr.Analog.AnalogData(:,fChan) = [];
    
end
if ~isempty(fLabel)
    Datastr.Force.ForceDataLabel = fLabel;
end
if isfield(Datastr,'Force')
    Datastr.Force = orderfields(Datastr.Force);
end

% EMG data
if ~isempty(eChan)
    disp(eChan)
    disp(fChan)
    Datastr.EMG.EMGData = Datastr.Analog.AnalogData(:,eChan);
    Datastr.EMG.EMGFrameRate = Datastr.Analog.AnalogFrameRate;
    
    %Datastr.Analog.AnalogData(:,eChan) = [];
    
end
if ~isempty(eLabel)
    Datastr.EMG.EMGDataLabel = eLabel;
end
if isfield(Datastr,'EMG')
    Datastr.EMG = orderfields(Datastr.EMG);
end

% Other data 
% This is currently somewhat obsolete, you might as well leave it in the analog channels
if ~isempty(oChan)
    
    Datastr.Other.OtherData = Datastr.Analog.AnalogData(:,oChan);
    Datastr.Other.OtherFrameRate = Datastr.Analog.AnalogFrameRate;
    
    %Datastr.Analog.AnalogData(:,oChan) = [];
    
end
if ~isempty(oLabel)
    Datastr.Other.OtherDataLabel = oLabel;
end
if isfield(Datastr,'Other')
    Datastr.Other = orderfields(Datastr.Other);
end

end