function [Datastr] = M5_getClusterAngles(Datastr,statTrial,outMode)
% gBMPDynUI statTrial = 1; outMode = 1;
% 
% INPUT)
% - Datastr: data structure
% - statTrial: integer or string, with the unique static trial identifier
% - outMode: string, either 'E','R', or 'ER'. When E is specified, output
% will be in Euler angles. When R is specified, output will be a 3x3
% rotation matrix. If ER is specified, both will be computed.
% 
% OUTPUT)
% Datastr: data structure, with added fields
%           .Marker.MarkerDataStatic
%           .Marker.ClusterAngleE
%           .Marker.ClusterAngleR
% 
% The latter two depending on which are specified with outMode.
% 
% NOTES)
% Requires you to have the static imported into a Datastr structure and
% stored into the subjroot folder
% 
% Requires soder.m


%% Check input

if isnumeric(statTrial)
    statTrial = num2str(statTrial);
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

%% Get median static marker data

Datastr = getStatic(Datastr,Datastrstat.Datastr);

%% Get cluster rotations

Datastr = getClusterAnglesParf(Datastr,outMode);

end