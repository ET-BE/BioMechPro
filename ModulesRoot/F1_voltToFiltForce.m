function [Datastr]  = F1_voltToFiltForce(Datastr,zeroFolder,zeroTrial,fZeroChan,MGRF,fType,fOrd,fCut)
% gBMPDynUI zeroFolder=1; zeroTrial=1; fZeroChan=1; MGRF=15; fType=1; fOrd=1; fCut=1;
%
% INPUT)
% - Datastr: structure, containing at least the fields:
% .Info.SubjRoot
% .Force.ForceData
%
% -zeroFolder: string, specifying the folder name relative to the root
% folder, that contains the zero measurement for the force data.
%
% -zeroFile: string, specifying a unique identifier of the zero trail, at
% the end of the file.
%
% - fZeroChan: scalar or vector, containing the channel numbers of the
% force data in the zero measurement, to be subtracted from the force data
% in the data structure.
%
% - MGRF: matrix, force plate ground reaction force calibration matrix.
%
% OUTPUT)
%
% NOTES)
%


%% Some checks

if isnumeric(zeroTrial)
    zeroTrial = num2str(zeroTrial);
end

if ~isfield(Datastr,'Force')
    warning('No field Force. Skipping.');
    return;
end

%% Get info

subjrootfolder = Datastr.Info.SubjRoot;
fs_frc = Datastr.Force.ForceFrameRate;

%% Find zero measurement file

if ~(isempty(zeroFolder) || isempty(zeroTrial))   % only do this when zero file is specified (for qualysis not necessarily needed)
    
    % Get all .C3D filenames in marker folder
    zerofiles = dir([subjrootfolder '\' zeroFolder '\*.mat']);
    zerofilenames = cell(1,length(zerofiles));
    zerofilenamescmp = cell(1,length(zerofiles));
    for ifile = 1:length(zerofiles)
        zerofilenames{ifile} = zerofiles(ifile).name;
        
        % Reverse string and remove '.mat', for comparison with strncmpi
        zerofilenamescmp{ifile} = zerofiles(ifile).name(end-4:-1:1);
    end
    clear zerofiles;
    
    % Check if storage filenames specified by zeroTrial field exist.
    % If not, skip the thing
    match = strncmpi( zeroTrial(end:-1:1),zerofilenamescmp,length(zeroTrial) );
    
    if sum(match) == 1
        zerofile = zerofilenames{match};
    else
        warning(['Skipping zero file ' zeroTrial ': not found or multiple found.']);
        return;
    end
    
    % Load XPCzero
    load([subjrootfolder '\' zeroFolder '\' zerofile])
end


%% Get filtered force data

% Get forces
fData = Datastr.Force.ForceData;
%fDataQTM = Datastr.ForceQTM.ForceData;

if isempty(MGRF)
    MGRFfile = ([subjrootfolder '\MGRF.mat']);
    load(MGRFfile);
end
if ~(isempty(zeroFolder) || isempty(zeroTrial))
    fData = ( fData - repmat(XPCzero(fZeroChan),[size(fData,1),1]) ) * MGRF;
else  % do not subtract zero file
    fData =  fData  * MGRF';
end

if ~(isempty(fOrd)||isempty(fCut)||isempty(fType))
    
    % Create filter
    [b,a] = butter(fOrd,2*fCut/fs_frc,fType);
    
    % Removing NaN values at end of fData
    ForceNaN = isnan(fData);
    ForceNonNaN = ~ForceNaN;
    for iforce = 1:12
        fDataNew(:,iforce) = fData(ForceNonNaN(:,iforce),iforce);
        %fDataQTMNew(:,iforce) = fDataQTM(ForceNonNaN(:,iforce),iforce);
    end
    fData = fDataNew;
    %fDataQTM = fDataQTMNew;
    % Filter force data
    fData = filtfilt(b,a,fData);
    %fDataQTM = filtfilt(b,a,fDataQTM);
end

% Rotation of the force data in the x y z configuration of OpenSim 
% (x -> AP, y -> up, z -> ML)
RotMat = [0 -1 0; 0 0 1; 1 0 0];
Ranges = [1:3; 4:6; 7:9; 10:12];
fData = fData';
%fDataQTM = fDataQTM';

for iPart = [1 2 3 4]
    iRange = Ranges(iPart,:);
    fData(iRange,:) = RotMat * fData(iRange,:);
    %fDataQTM(iRange,:) = RotMat * fDataQTM(iRange,:);
end
fData = fData';
%fDataQTM = fDataQTM';

% Store
Datastr.Force.ForceData = fData;
%Datastr.ForceQTM.ForceData=fDataQTM;

end