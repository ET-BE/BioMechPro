function [Datastr] = M3_getProbe(Datastr,probefolder,probeMarkerLocal,idxProbe,clusterCell,probeCell)
% gBMPDynUI probeFolder=1; pMarkLocal=5; idxProbe=1; clusterCell=3; probeCell=10;
%
% The same transformation and filter settings are applied as on the marker
% data from the trial, if this info is available.
%
% INPUT)
% - Datastr : the data structure, with no specific required fields
%
% - probefolder : string, specifying the folder containing the probe trials,
% relative to the subject rootfolder. Example:
% 'VZ\Probe'
%
% - probeMarkerLocal : matrix, containing the distances in meters of the
% probe LEDs relative to the probe tip. Example (UTWENTE probe):
% [ 120+28.65      45     10.5 ; ...
%   120+28.65     -45     10.5 ; ...
%   195.35+28.65   75     10.5 ; ...
%   195.35+28.65  -75     10.5 ; ...
%   28.65           0     10.5] ./1000;
% Which contains the distances from the tip to the lower left, lower right,
% upper left, upper right, and tip LEDs, in meters.
%
% - idxProbe : vector, with indices specifying the probe LEDs in the marker
% data, in the same order as specified in probeMarkerLocal.
% Example :
% [20 21 23 19 24];
% Where 20=lower left, 21=lower right 23=upper left 19=upper right, 24=tip
%
% - clusterCell : cell, containing the unique names given to the marker
% clusters that were used. Example :
% {'LFOOT','RFOOT','LLOWER','RLOWER','LUPPER','RUPPER','PELVIS','TRUNK','HEAD'};
%
% - probeCell : cell array, specifying which probe trials belong to which
% cluster in the clusterCell. probeCell contains the filenames of the probe
% trials. Each row in probeCell corresponds with every entry in clusterCell
% Use '' to ensure the array has the same number of columns in every row.
% Example:
% {'CAL','CM1L','CM5L','';...
%  'CAR','CM1R','CM5R','';...
%  'MLL','MML','CFL','';...
%  'MLR','MMR','CFR','';...
%  'ELL','EML','TML','';...
%  'ELR','EMR','TMR','';...
%  'SIASL','SIPSL','SIASR','SIPSR';...
%  'C7','PX','IJ','';...
%  'OCC','HV','SEL',''};
% In which the rows belong to the LFOOT, RFOOT, LOWER, RLOWER, LUPPER,
% RUPPER, PELVIS, TRUNK and HEAD marker clusters respectively.
%
% OUTPUT)
% - Datast : the data structure, with added fields:
% .Marker.ProbeDataLocal.ProbePos.PosName
% .Marker.ProbeDataLocal.ProbePos.ClusterLabel
% .Marker.ProbeDataLocal.ProbePos.UsedLEDs
% .Marker.ProbeDataLocal.ProbePos.ClusterMarkerLocal
%
% NOTES)
%

%% Get info

mTrans = Datastr.Marker.MarkerDataTransform;

if isfield(Datastr.Marker,'MarkerDataFilter')
    filterCoeff = Datastr.Marker.MarkerDataFilter;
else
    filterCoeff = [];
end

%% Local Probe information

% Get root folder
rootfolder = Datastr.Info.SubjRoot;

n = 0; % Second loop variable
for i = 1:length(clusterCell)
    
    % Select the probe positions belonging to each cluster
    probePerCluster = probeCell(i,:);
    
    for j = 1:length(probePerCluster)
        if ~isempty(probePerCluster{j})
            n = n + 1;
            
            % Load probe file
            FullFileName = [rootfolder '\' probefolder '\' probePerCluster{j} '.C3D'];
            [C3DRaw] = readC3D(FullFileName);
            
            % Do basic marker processing
            C3DProc = procC3D(C3DRaw,mTrans); % Mind that this interpolates data
            
            % Apply zero phase filter, if any
            if ~isempty(filterCoeff)
                C3DProc.Marker.MarkerData = filtfilt(filterCoeff{1},filterCoeff{2},C3DProc.Marker.MarkerData);
            end
            
            % Operations based on 'keepInput', that also need to be applied to the probe data.
            if isfield(Datastr.Info,'zpbutterMarker2')
                fType = Datastr.Info.zpbutterMarker2;
                fOrd = Datastr.Info.zpbutterMarker3;
                fCut = Datastr.Info.zpbutterMarker4;
<<<<<<< HEAD
                C3DProc = M2_zpbutterMarker(C3DProc,fType,fOrd,fCut,[]);
=======
                C3DProc.Info = Datastr.Info;
                C3DProc = M2_zpbutterMarker(C3DProc,fType,fOrd,fCut,[]); 
>>>>>>> origin/extraModules
            end
            if isfield(Datastr.Info,'editLabel2')
                labelField = Datastr.Info.editLabel2;
                cellStrOld = Datastr.Info.editLabel3;
                cellStrNew = Datastr.Info.editLabel4;
                C3DProc = Gen_editLabel(C3DProc,labelField,cellStrOld,cellStrNew,[]);
            end
            if isfield(Datastr.Info,'editMarker2')
                dataField = Datastr.Info.editMarker2;
                dim = Datastr.Info.editMarker3;
                markIdxOld = Datastr.Info.editMarker4;
                markIdxNew = Datastr.Info.editMarker5;
                
                % NOTE: assumption is here that all probe markers follow up
                for idx = 1:length(markIdxOld)
                    if any(markIdxOld(idx) >= idxProbe)
                        markIdxOld(idx) = markIdxOld(idx) + numel(idxProbe);
                    end
                    if any(markIdxNew(idx) >= idxProbe)
                        markIdxNew(idx) = markIdxNew(idx) + numel(idxProbe);
                    end
                end
                
                C3DProc = Gen_editMarker(C3DProc,dataField,dim,markIdxOld,markIdxNew,[]);
            end
            
            % Get local information of the probe position w.r.t. the marker frame
<<<<<<< HEAD
            %             disp([probePerCluster{j} ' ' clusterCell{i}]) % For debugging / when encountered error below
=======
             disp([probePerCluster{j} ' ' clusterCell{i}]) % For debugging / when encountered error below
>>>>>>> origin/extraModules
            [clusterMarkerLocal,usedLEDs] = getProbeLocal(C3DProc,clusterCell{i},idxProbe,probeMarkerLocal);
            
            % Store all probe information for later use
            Body.ProbePos{n}.PosName = probePerCluster{j};
            Body.ProbePos{n}.ClusterLabel = clusterCell{i};
            Body.ProbePos{n}.UsedLEDs = usedLEDs;
            Body.ProbePos{n}.ClusterMarkerLocal = clusterMarkerLocal;
            % (you could opt to store this in the Datastr structure as well)
            
            
        end
    end
end

%% Reconstruct global probes

Datastr = getProbeGlobalParf(Datastr,Body);

end