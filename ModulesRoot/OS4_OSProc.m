function [Datastr] = OS4_OSProc(Datastr,statTrial)
% gBMPDynUI statTrial=1;
% 
% Do some processing on the OS data once it has been imported into Matlab
% (Remove mean static angles, optional: correction for strange joint angles)

% Some checks
if isnumeric(statTrial)
    statTrial = num2str(statTrial);
end

% Do a correction for weird joint angle shifts
% Datastr = correctIKAngle(Datastr);

% Find and load the static
root = Datastr.Info.SubjRoot;
backslsh = strfind(root,'\');
fullfilename = [root '\' root(backslsh(end)+1:end) '_' statTrial '.mat'];
foo = load(fullfilename);

DatastrStat = foo.Datastr;
clear foo; 

% Remove mean static angles
datastat = DatastrStat.Marker.JointAngData;
data = Datastr.Marker.JointAngData;

Datastr.Marker.JointAngData = data - repmat(median(datastat), [size(data,1) 1 1] );


end