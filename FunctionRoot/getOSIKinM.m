function [C3Ddata] = getOSIKinM(C3Ddata,subjosikmot,varargin)
%% getOSIDinM
% 
% Import OpenSim inverse dynamics joint torque data from a .sto file into 
% the corresponding C3DProc data structure.
% 
% INPUT)
% C3Ddroc : C3Ddata structure
% subjosikmot : string - filename of the .mot file containing the IK joint angles.
% permvec : permutation vector to put data in C3DProc coordinates
% In OS the directions are defined as: x: walking direction, z: to the right, y: upward)
% Default is [1 2 3] (no permutation)
% 
% OUTPUT)
% C3Ddata : C3Ddata strucure with added joint angles
% 
% NOTES)
% This file might require expansion for other joint angles

% Update history
% 25-09-2015 : Mark Vlutters : File creation


%% Check input

if nargin < 3  || isempty(varargin{1})
    permvec = [1 2 3];
else
    permvec = varargin{1};
end


%% Load data

% Read file header and data header
fid = fopen(subjosikmot,'r');
if fid == -1
    error('getOSIKinM:fileopen',['Failed to open ' subjosikmot]);
end

tline = '';
ntimeout = 0;
while isempty(strfind(tline,'time')) && ntimeout < 20
    
    tline = fgetl(fid);
    
    if ~isempty(strfind(tline,'nRows='))
        nRows = str2double( tline(strfind(tline,'nRows=') + 6 : end) );
    elseif ~isempty(strfind(tline,'nColumns='))
        nCols = str2double( tline(strfind(tline,'nColumns=') + 9 : end) );
    elseif ~isempty(strfind(tline,'inDegrees='))
        isDeg = tline(strfind(tline,'inDegrees=') + 10 : end);
    end
    
    ntimeout = ntimeout + 1;
end
fclose(fid);

if ntimeout == 20
    error('getOSIKinM:DataHeader','Unable to find data column header in file.');
else
    angData = dlmread(subjosikmot,'\t',[ntimeout,0,ntimeout+nRows-1,nCols-1]);
    % NOTE: assumed here the data is directly below the column header
end


% %% Sort data based on header
% % 
% % NOTE: It is assumed the column header is delimited with tabs
% % NOTE2: Even though the torques are not in global coordinates, they are
% % sorted in a "global as possible" way. For example, the knee
% % flexion-extension moment is assumed to be in the sagittal plane, about
% % the transversal axis, regardless of e.g. the upper leg rotation.
% % As a reminder, in OpenSim, x is the walking direction, y is upward and z
% % to the right
% 
headercell = strsplit(tline,'\t');

% % Column and layer index for data matrix
% stoidx = zeros(length(headercell),3);
% 
% % Create indices and labelcells for data storage based on header names
% for ihd = 1:length(headercell)
%     switch headercell{ihd}
%         case 'pelvis_tx'
%             stoidx(ihd,:) = [1 1 ihd];
%         case 'pelvis_ty'
%             stoidx(ihd,:) = [1 2 ihd];
%         case 'pelvis_tz'
%             stoidx(ihd,:) = [1 3 ihd];
%         case 'pelvis_list'
%             stoidx(ihd,:) = [2 1 ihd];
%         case 'pelvis_rotation'
%             stoidx(ihd,:) = [2 2 ihd];
%         case 'pelvis_tilt'
%             stoidx(ihd,:) = [2 3 ihd];
%           
%         case 'lumbar_bending'
%             stoidx(ihd,:) = [3 1 ihd];
%         case 'lumbar_rotation'
%             stoidx(ihd,:) = [3 2 ihd];
%         case 'lumbar_extension'
%             stoidx(ihd,:) = [3 3 ihd];
%         
%         case 'hip_adduction_l'
%             stoidx(ihd,:) = [4 1 ihd];
%         case 'hip_rotation_l'
%             stoidx(ihd,:) = [4 2 ihd];
%         case 'hip_flexion_l'
%             stoidx(ihd,:) = [4 3 ihd];
%         case 'knee_angle_l'
%             stoidx(ihd,:) = [5 3 ihd];
%         case 'subtalar_angle_l'
%             stoidx(ihd,:) = [6 1 ihd];
%         case 'ankle_angle_l'
%             stoidx(ihd,:) = [6 3 ihd];
%         
%         case 'hip_adduction_r'
%             stoidx(ihd,:) = [7 1 ihd];
%         case 'hip_rotation_r'
%             stoidx(ihd,:) = [7 2 ihd];
%         case 'hip_flexion_r'
%             stoidx(ihd,:) = [7 3 ihd];
%         case 'knee_angle_r'
%             stoidx(ihd,:) = [8 3 ihd];
%         case 'subtalar_angle_r'
%             stoidx(ihd,:) = [9 1 ihd];
%         case 'ankle_angle_r'
%             stoidx(ihd,:) = [9 3 ihd];
%             
%         otherwise % time / metatarsophalangeal joint
%             % Do nothing
%     end
% end
% 
% % Remove unassigned indices from stoidx (e.g. if not all present in .sto file)
% stoidx(sum(stoidx,2) == 0,:) = [];
% 
% % Create corresponding labelcell
% % NOTE: you should expand this cell if you want to add other joints.
% % Mind the order of the labels. They should correspond with the number order in the first column of stoidx.
% labelcellAll = {'PELVIST','PELVISR','LJC','HJCL','KJCL','AJCL','HJCR','KJCR','AJCR'};
% convertMask = [false true true true true true true true true]; % For conversion angle -> radian
% 
% iCols = unique(stoidx(:,1));
% jointAngDataLabel = labelcellAll(iCols);
% convertMask = convertMask(iCols);

%% Assign data

% Make data equally long as the rest of the marker data using the sync indices
% msiz = size(C3Ddata.Marker.MarkerData,1);
% if isfield(C3Ddata.Marker,'MarkerSyncIdx')
%     syncIdx = C3Ddata.Marker.MarkerSyncIdx;
% else
%     syncIdx = [1 msiz];
% end
% 
% % Pre-alloc new data container
% jointAngData = zeros(msiz,numel(iCols),3);

% Put the data in a new matrix
% for ihd = 1:size(stoidx,1)
%     jointAngData(syncIdx(1):syncIdx(1)+size(angData,1)-1,stoidx(ihd,1),stoidx(ihd,2)) = angData(:,stoidx(ihd,3));
% end
% 
% % Permutate the data dimensions
% jointAngData = jointAngData(:,:,permvec);
% 
% % Convert to radians if in degrees
% if strcmpi(isDeg,'yes')
%      jointAngData(:,convertMask,:) = jointAngData(:,convertMask,:) .* pi/180;
% end

% Add data to C3DProc structure
C3Ddata.Resample.Sych.IKAngData = angData;
C3Ddata.Resample.Sych.IKAngDataLabel = headercell;

end