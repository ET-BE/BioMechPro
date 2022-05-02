function [Datastr] = A1_sortAnalogX(Datastr,xanaFolder,fChan,fLabel,eChan,eLabel,oChan,oLabel)
% gBMPDynUI xanaFolder=1; fChan=1; fLabel=1; eChan=1; eLabel=1; oChan=1; oLabel=1;
% 
% Import and sort external analog data

%% Find XPC file

if ~isempty(fChan) || ~isempty(eChan) || ~isempty(oChan) % In case you only want to swap labels

    rootfolder = Datastr.Info.SubjRoot;
    subjtrial = Datastr.Info.Trial;

    % Get all .mat filenames in XPC folder
    xpcfiles = dir([rootfolder '\' xanaFolder '\*.mat']);
    xpcfilenames = cell(1,length(xpcfiles));
    xpcfilenamescmp = cell(1,length(xpcfiles));
    for ifile = 1:length(xpcfiles)
        xpcfilenames{ifile} = xpcfiles(ifile).name;

        % Reverse string and remove '.mat', for comparison with strncmpi
        xpcfilenamescmp{ifile} = xpcfiles(ifile).name(end-4:-1:1); 
    end
    
    % Check if storage filenames specified by Trials field exist. 
    % If not, skip the thing
    match = strncmpi( subjtrial(end:-1:1),xpcfilenamescmp,length(subjtrial) );
    
    if sum(match)==1
        subjxpcfile = xpcfilenames{match};
    else
        warning('A1_sortAnalogX:nofile',['Could not find unique XPC file ' subjtrial '. Skipping.']);
        return;
    end

    fullFileName = [rootfolder '\' xanaFolder '\' subjxpcfile]; 

    load(fullFileName);

end

%% Get XPC data

% Force data
if ~isempty(fChan)
    Datastr.Force.ForceData = XPCdata.data(:,fChan);
    Datastr.Force.ForceFrameRate = round( 1 ./ mean(diff(XPCdata.data(:,end))) ); % Assumed last channel is XPC time vector
end
if ~isempty(fLabel)
    Datastr.Force.ForceDataLabel = fLabel;
end
if isfield(Datastr,'Force')
    Datastr.Force = orderfields(Datastr.Force);
end

% EMG data
if ~isempty(eChan)
    Datastr.EMG.EMGData = XPCdata.data(:,eChan);
    Datastr.EMG.EMGFrameRate = round( 1 ./ mean(diff(XPCdata.data(:,end))) ); % Assumed last channel is XPC time vector
end
if ~isempty(eLabel)
    Datastr.EMG.EMGDataLabel = eLabel;
end
if isfield(Datastr,'EMG')
    Datastr.EMG = orderfields(Datastr.EMG);
end

% Other data
if ~isempty(oChan)
    Datastr.Other.OtherData = XPCdata.data(:,oChan);
    Datastr.Other.OtherFrameRate = round( 1 ./ mean(diff(XPCdata.data(:,end))) ); % Assumed last channel is XPC time vector
end
if ~isempty(oLabel)
    Datastr.Other.OtherDataLabel = oLabel;
end
if isfield(Datastr,'Other')
    Datastr.Other = orderfields(Datastr.Other);
end

end