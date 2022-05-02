function [DatastrReset] = R1_strucReset(Datastr)
% Reset data structure to Marker and Probe data only
% Note that filter operations are not undone.

%% Do stuff

DatastrReset.Marker.MarkerData = Datastr.Marker.MarkerData;
DatastrReset.Marker.MarkerDataLabel = Datastr.Marker.MarkerDataLabel;
DatastrReset.Marker.MarkerFrameRate = Datastr.Marker.MarkerFrameRate;
DatastrReset.Marker.ProbedData = Datastr.Marker.ProbedData;
DatastrReset.Marker.ProbedDataLabel = Datastr.Marker.ProbedDataLabel;
DatastrReset.Marker.ProbedDataLabelC = Datastr.Marker.ProbedDataLabelC;

DatastrReset.Analog.AnalogData = Datastr.Analog.AnalogData;
DatastrReset.Analog.AnalogDataLabel = Datastr.Analog.AnalogDataLabel;
DatastrReset.Analog.AnalogFrameRate = Datastr.Analog.AnalogFrameRate;


end