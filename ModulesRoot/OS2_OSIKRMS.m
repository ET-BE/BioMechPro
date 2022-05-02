function [Datastr] = OS2_OSIKRMS(Datastr,osFolder)
% gBMPDynUI osFolder=1;
% 
% Calculate the marker RMS error of the IK, per subject over all samples
% per trial. Output will be appended to a file which displays the overall
% results.
%
% INPUT)
% - Datastr: with no specific fields required
% 
% - osFolder: string, specifying the folder in which the OpenSim files will
% be stored, relative to the subject root folder. Example:
% 'OS'
% 
% OUTPUT)
% - No direct output
% - A file named OSIKRMS.m will be written to the current path, containing
% all RMS values from all trials.
% 
% NOTES)


%% Get info and paths

rootfolder = Datastr.Info.SubjRoot;
bckslsh = strfind(rootfolder,'\');
if isempty(bckslsh)
    bckslsh = 0;
end

trial = Datastr.Info.Trial;
savename = [rootfolder(bckslsh(end)+1:end) trial];

ikLogPath = [rootfolder '\' osFolder '\' savename 'IKout.log'];
writePath = [pwd '\OSIKRMS.m'];

%% Calculate error values

% Get values
[rmsm] = getOSIKRMS(ikLogPath,'median');
% [rmsm] = getOSIKRMS(ikLogPath,'mean');

% Write to file
if ~exist(writePath,'file')
    fid = fopen(writePath,'w+');
else
    fid = fopen(writePath,'a+');
end
fprintf(fid,'%f \n',rmsm);
fclose(fid);

% Set Datastr to empty to prevent saving
Datastr = [];

end
