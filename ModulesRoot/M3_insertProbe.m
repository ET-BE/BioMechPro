function [Datastr] = M3_insertProbe(Datastr,statTrial,probeNames,clusterNames,insertIdx)
% gBMPDynUI statTrial=1; probeNames=2; clusterNames=2; insertIdx=1;
% 
% Manually insert probe points.
% 
% INPUT)
% 
% OUTPUT)
% 
% NOTES)
% This is an ugly piece of code. Try not to use this unless you seriously
% messed up and forgot probe points.


%% Check input

if isnumeric(statTrial)
    statTrial = num2str(statTrial);
end

if ischar(probeNames)
    probeNames = {probeNames};
end

if ischar(clusterNames)
    clusterNames = {clusterNames};
end

%% Find and load static file

subjrootfolder = Datastr.Info.SubjRoot;

% Get all .mat filenames in rootfolder
matfiles = dir([subjrootfolder '\*.mat']);
matfilenames = cell(1,length(matfiles));
matfilenamescmp = cell(1,length(matfiles));
for ifile = 1:length(matfiles)
    matfilenames{ifile} = matfiles(ifile).name;

    % Reverse string and remove '.mat', for comparison with strncmpi
    matfilenamescmp{ifile} = matfiles(ifile).name(end-4:-1:1); 
end
clear matfiles;

% Check if storage filenames specified by Trials field exist. 
% If not, skip the thing
match = strncmpi( statTrial(end:-1:1),matfilenamescmp,length(statTrial) );
if sum(match) == 1
    subjmatfile = matfilenames{match};
else
    warning(['Skipping static trial' statTrial ': not found or multiple found.']);
    return;
end

fullfilename = [subjrootfolder '\' subjmatfile];

Datastrstat = load(fullfilename);


%% Get probe information

tmpname = Datastr.Info.SubjRoot(end-3:end);

if ~exist(['M3_insertProbe_Body_TMP_' tmpname '.mat'],'file')

    for iprb = 1:length(probeNames)

        [clusterMarkerLocal,usedLEDs] = putProbe(Datastrstat.Datastr,clusterNames{iprb});

         % Store all probe information for later use
        Body.ProbePos{iprb}.PosName = probeNames{iprb};
        Body.ProbePos{iprb}.ClusterLabel = clusterNames{iprb};
        Body.ProbePos{iprb}.UsedLEDs = usedLEDs;
        Body.ProbePos{iprb}.ClusterMarkerLocal = clusterMarkerLocal;
        % (you could opt to store this in the Datastr structure as well)

    end

    save(['M3_insertProbe_Body_TMP_' tmpname],'Body')
    
else
    
    load(['M3_insertProbe_Body_TMP_' tmpname '.mat']);
    
end

%% Reconstruct global probes

foo = Datastr;
foo.Marker = rmfield(Datastr.Marker,{'ProbedData','ProbedDataLabel','ProbedDataLabelC'});
foo = getProbeGlobalParf(foo,Body);

%% Merge inserted probes with actual data structure
% Using placement info provided by user

% Allocate
ProbedDataNew = zeros(size(Datastr.Marker.ProbedData,1),size(Datastr.Marker.ProbedData,2)+length(insertIdx),3);
ProbedDataLabelNew = cell(1,length(Datastr.Marker.ProbedDataLabel)+length(insertIdx));
ProbedDataLabelCNew = cell(1,length(Datastr.Marker.ProbedDataLabel)+length(insertIdx));

% Assign new
ProbedDataNew(:,insertIdx,:) = foo.Marker.ProbedData;
ProbedDataLabelNew(insertIdx) = probeNames;
ProbedDataLabelCNew(insertIdx) = clusterNames;

% Assign old
n = 1; % Additional counter
for idx = 1:size(ProbedDataNew,2)

    if sum(sum(ProbedDataNew(:,idx,:),1),3) == 0
        ProbedDataNew(:,idx,:) = Datastr.Marker.ProbedData(:,n,:);
        ProbedDataLabelNew{idx} = Datastr.Marker.ProbedDataLabel{n};
        ProbedDataLabelCNew{idx} = Datastr.Marker.ProbedDataLabelC{n};
        n = n + 1;
    end
    
end

%% Assign to final structure

Datastr.Marker.ProbedData = ProbedDataNew;
Datastr.Marker.ProbedDataLabel = ProbedDataLabelNew;
Datastr.Marker.ProbedDataLabelC = ProbedDataLabelCNew;

end