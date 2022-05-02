function [rmsm,rmsstd,rmsval] = getOSIKRMS(subjosiklog,mode)
%% getOSIKRMS
% Get mean and standard deviation of the OpenSim inverse kinematics RMS 
% error read from the IKout.log file.
% 
% INPUT)
% subjosiklog : string - filename of the .log file containing the output
% log of the IK, including .log extension.
% 
% mode : string - either 'mean' or 'median'. Determines the output mode for
% rmsm.
% 
% OUTPUT)
% rmsm : mean RMS over all samples
% rmsstd : standard deviation of the RMS of all samples
% 
% NOTES)
% 

% Update history)
% 02-10-2015 - Mark Vlutters - File creation


%% Read file and get data

% Open file
fid = fopen(subjosiklog,'r');
if fid == -1
    error('getIKRMS:cannotopen',['Failed to open ' subjosiklog]);
end

% Read file
tline = '';
rmsval = zeros(60000,1); % Pre-alloc for 10 min at 100 Hz
loopc = 0;
while ~isnumeric(tline)
    tline = fgetl(fid);
    linidx = strfind(tline,'RMS=');
    if ~isempty(linidx)
        loopc = loopc + 1;
        rmsval(loopc) = str2double(strtok(tline(linidx+4:end),','));
    end
end

% Close file
fclose(fid);

% Clean up rmsval
rmsval(rmsval == 0) = [];

% Calculate m and std
if strcmpi(mode,'median')
    rmsm = median(rmsval); % Deals better with outlier samples with high error
else
    rmsm = mean(rmsval);
end
rmsstd = std(rmsval,[],1);

end