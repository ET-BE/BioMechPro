function getCEINMStrials(osmaFilePath,osidstoFilePath,trialsGenSetPath,emgFilePath)


xmlSet = xmlread([ trialsGenSetPath ] ); % 
% Modify settings in xml
xmlSet.getElementsByTagName('muscleTendonLengthFile').item(0).setTextContent([osmaFilePath '\' '_MuscleAnalysis_Length.sto']); %
xmlSet.getElementsByTagName('excitationsFile').item(0).setTextContent(emgFilePath); %
xmlSet.getElementsByTagName('momentArmsFile').item(0).setTextContent([osmaFilePath '\' '_MuscleAnalysis_MomentArm_knee_angle_r.sto']); %
xmlSet.getElementsByTagName('momentArmsFile').item(1).setTextContent([osmaFilePath '\' '_MuscleAnalysis_MomentArm_ankle_angle_r.sto']); %
xmlSet.getElementsByTagName('externalTorquesFile').item(0).setTextContent(osidstoFilePath); %

% Get info from paths
out1=regexp(osmaFilePath,'\','split');
mocapfolder=out1{end-1};
subject=out1{end};
subject=subject(1:end-2);
filename=[subject '.xml'];

out2=regexp(trialsGenSetPath,'\','split');
pathCEINMS=out2(1:end-1);
pathCEINMS=fullfile(pathCEINMS{:});

%Save xml
setFile = [pathCEINMS '\' mocapfolder '\' filename] ;
xmlwrite(setFile,xmlSet);