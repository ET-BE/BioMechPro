function [] = getOSIK(osinstallpath,subjosmod,subjosikset,subjostrc)
%% getOSIK
% Do inverse kinematics (IK) using OpenSim
% 
% INPUT)
% osinstallpath : string - OpenSim installation path (e.g. C:\Program Files\OpenSim\OpenSim 3.3)
% subjosmod : string - .osim file name containing the subject specific
%           OpenSim model used for the inverse kinematics.
% subjosikset : string - generic .xml setup file name for OpenSim inverse kinematics, with .xml extension
% subjostrc : string - .trc file name containing the marker data on which
%           the inverse kinematics must be applied, with .trc extension.
% 
% For the last 3 inputs, if no full path is supplied the file is assumed 
% present in the current folder.
% 
% OUTPUT)
% No direct output is supplied.
% The code generates 3 or 4 output files:
% - A trial specific .xml IK settings file
% This file is saved in the same folder as the generic .xml file
% - A .mot file containing the IK output
% This file will contain 1 sample less than the .trc file (last sample)
% This file is saved in the same folder as the .trc marker file
% - Two .log files
% out.log will contain the console output, err.log is an error log.
% The error log will be deleted if it is empty.
% These files are moved to the same folder as the .trc marker file
% 
% NOTES)
% Before being able to run this code you need:
% .osim model containing an OpenSim model used in the IK
% .trc files containing marker data on which the IK must be applied
% .xml file containing generic IK settings
%
% The .osim model is preferably a scaled, subject specific model. It is
% good practice to take note of the scaling output. (RMS < 0.01).
% 
% It is assumed that the .trc file contains the time samples in the 2nd
% column of the file, starting at the 7th row.
% 
% The generic .xml settings file must be created in OpenSim by creating one
% using the save button in the Inverse Kinematics tool. 
% 
% This code was tested with OpenSim 3.3
% 
% Tip: for more information on xml DOM editing see: 
% http://www.w3schools.com/xml/dom_node.asp

% Update history)
% 22-09-2015 : Mark Vlutters : file creation

%% Read and prepare settings

% Read sample numbers from .trc file
trcData = dlmread(subjostrc,'\t',6,1); 
% NOTE: dlmread uses zero indexing (7th row, 2nd column)
% NOTE2: this is inefficient, but since the number of rows in the
% .trc file is not available here beforehand, you have to read all data 
% from the upper left corner.

% Create .trc output file path
% If full path is supplied, take last part for name inside .trc file
if ~isempty(strfind(subjostrc,'\'))
    foo = strfind(subjostrc,'\');
    trcfilename = regexprep(subjostrc(foo(end)+1:end),'.trc','');
    trcfilepath = subjostrc(1:foo(end-1));
else
    trcfilename = regexprep(trialname,'.trc');
    trcfilepath = pwd;
end

% Read generic settings file
xmlSet = xmlread(subjosikset);

% Create settings output file path
% If full path is supplied, take last part for name inside trc file
if ~isempty(strfind(subjosikset,'\'))
    foo = strfind(subjosikset,'\');
%     setfilename = subjosikset(foo(end)+1:end);
    setfilepath = subjosikset(1:foo(end));
else
%     setfilename = trialname;
    setfilepath = pwd;
end

% Modify settings file to a trial specific one
% NOTE: For good practice: metadata (data about data) should be stored as attributes, and the data itself should be stored as elements.
% NOTE2: rather than setTextContent you could also use .getFirstChild.setNodeValue, but this only works if the node is already a text node.
% Beware that setTextContent removes all child nodes from the node to which it is applied, and replaces it with the supplied text. 
% (Although not relevant here as model_file, marker_file, output_motion_file and time_range have no other content).
% NOTE3: the assumption here is that these nodes already exist in the .xml file.
% Other verions of OpenSim might use different node names (e.g. in older
% versions of OpenSim the node InverseKinematicsTool was called IKTool).
if xmlSet.getElementsByTagName('InverseKinematicsTool').item(0).hasAttributes
    xmlSet.getElementsByTagName('InverseKinematicsTool').item(0).getAttributes.item(0).setValue(subjostrc); 
else
    xmlSet.getElementsByTagName('InverseKinematicsTool').item(0).setAttribute('name',subjostrc)
end
xmlSet.getElementsByTagName('model_file').item(0).setTextContent(subjosmod); % Input model file
xmlSet.getElementsByTagName('marker_file').item(0).setTextContent(subjostrc);  % Input marker file 
xmlSet.getElementsByTagName('output_motion_file').item(0).setTextContent([trcfilepath 'DataFiles\' trcfilename 'IK.mot']); % Output file
xmlSet.getElementsByTagName('time_range').item(0).setTextContent([num2str(trcData(1,1)) ' ' num2str(trcData(end,1))]); % Output range

clear trcData;

% Write modified settings .xml file
setFile = [setfilepath trcfilename 'IKset.xml'];
xmlwrite(setFile,xmlSet);

%% Do inverse kinematics

% Send system command to perform IK using specific settings
disp(['Starting IK for ' trcfilename]);
system(['"' osinstallpath '\bin\ik.exe" -Setup ' setFile ' > nul']);
% NOTE: the > nul suppresses the window output

% Rename generic IK output log and move it to the right folder

% Check if folder exist, if not create new
if ~ exist([trcfilepath '\Logs\'],'dir')
    mkdir([trcfilepath '\Logs\']);
end

if exist('out.log','file')
    movefile('out.log',[trcfilepath '\Logs\' trcfilename 'IKout.log']); % Output log
end
if exist('err.log','file')
    
    % Check if error log is empty
    fid = fopen('err.log','r');
    if fid ~= -1
        foo = fread(fid);
        fclose(fid);
    else
        foo = true;
    end
    
    % Remove error log if empty, otherwise move and rename
    if isempty(foo)
        delete('err.log');
    else
        warning('getOSIK:errors',['Please see ' trcfilename 'IKerr.log']);
        movefile('err.log',[trcfilepath '\Logs\' trcfilename 'IKerr.log']); % Error log
    end
end

disp(['Finished IK for ' trcfilename '. See .log file for details.']);

end