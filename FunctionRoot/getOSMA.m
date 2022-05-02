function [] = getOSMA(osinstallpath,subjosmod,subjosmaset,subjosikmot)
%% getOSMA
% Do muscle analysis (MA) using OpenSim
% 
% INPUT)
% osinstallpath : string - OpenSim installation path (e.g. C:\Program Files\OpenSim\OpenSim 3.3)
% subjosmod : string - name of the .osim file containing the subject specific
%             OpenSim model to use for the inverse dynamics.
% subjosmaset : string - generic .xml setup file name for OpenSim  muscle analysis, with .xml extension
% subjosikmot : string - name of the .mot file containing motion data, e.g. as
%               generated from the inverse kinematics.
% subjosxldmot : string - name of the .mot file containing external load
%               data, such as ground reaction forces and moments
% 
% OUTPUT)
% No direct output is supplied.
% The code generates an .sto file containing variables selected in MA xml.

% 
% NOTES)
% Before running this file, generic settings files should be created using
% OpenSim. The general MA settings file can be saved directly from the MA
% tool screen. All output files are stored in a MA folder in the folder where the subjosmaset file is
% located. The filename of the IK.mot file is reused for naming of output
% files. 
% 
% It is assumed that the name of the subjosikmot file ends on IK.mot, 
% as generated by getOSIK.m
% 
% It is assumed that the subjosikmot file contains the time values in
% column 1, starting at row 12. If this is not the case you should adjust
% the indices in dlmread below.
% 
% If no full path is supplied the file is assumed present in the current folder.
% 
% 
% This code was tested with OpenSim 3.3

%% Get some directories

% Extract output directory from ma.xml file 
if ~isempty(strfind(subjosmaset,'\'))
    foo = strfind(subjosmaset,'\');
    stofilepath = subjosmaset(1:foo(end));
else
    stofilepath = pwd;
end

% Extract general output file name from IK.mot file
if ~isempty(strfind(subjosikmot,'\'))
    foo = strfind(subjosikmot,'\');
    genfilename = regexprep(subjosikmot(foo(end)+1:end),'IK.mot','');
else
    genfilename = regexprep(subjosikmot,'IK.mot','');
end


%% Read and prepare general settings

% Open IK.mot and check data length
ikmotdata = dlmread(subjosikmot,'\t',11,0); % NOTE: zero indexing

% General MA Settings
xmlSet = xmlread(subjosmaset);
% Modify settings to a trial specific one
xmlSet.getElementsByTagName('results_directory').item(0).setTextContent([stofilepath genfilename 'MA']); % Output path
xmlSet.getElementsByTagName('model_file').item(0).setTextContent(subjosmod); % Model file
xmlSet.getElementsByTagName('initial_time').item(0).setTextContent(num2str(ikmotdata(1,1)));
xmlSet.getElementsByTagName('final_time').item(0).setTextContent(num2str(ikmotdata(end,1)));
xmlSet.getElementsByTagName('start_time').item(0).setTextContent(num2str(ikmotdata(1,1)));
xmlSet.getElementsByTagName('end_time').item(0).setTextContent(num2str(ikmotdata(end,1)));
xmlSet.getElementsByTagName('coordinates_file').item(0).setTextContent(subjosikmot); % IK.mot

% Write modified settings .xml file
setFile =  [stofilepath genfilename 'MA' '\' 'MAset.xml'];
xmlwrite(setFile,xmlSet);

%% Do inverse dynamics

% Send system command to perform ID using specific settings
disp(['Starting MA for ' genfilename]);
system(['"' osinstallpath '\bin\analyze.exe" -Setup ' setFile ' > nul']);
% NOTE: the > nul suppresses the window output

% Rename generic MA output log and move it to the right folder
if exist('out.log','file')
    movefile('out.log',[stofilepath genfilename 'MA\' 'MAout.log']); % Output log
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
        warning('getOSID:errors',['Please see ' genfilename 'MA\' 'MAerr.log']);
        movefile('err.log',[stofilepath genfilename 'MA\' 'MAerr.log']); % Error log
    end
end

disp(['Finished MA for ' genfilename '. See .log file for details.' ]);


end