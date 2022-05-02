function [C3Ddata] = getOSIDinM_Portable(C3Ddata,subjosidsto,varargin)
%% getOSIDinM
% 
% Import OpenSim inverse dynamics joint torque data from a .sto file into 
% the corresponding C3DProc data structure.
% 
% INPUT)
% C3Ddroc : C3Ddata structure
% subjosidsto : string - filename of the .sto file containing the joint torques.
% permvec : permutation vector to put data in C3DProc coordinates
% In OS the directions are defined as: x: walking direction, z: to the right, y: upward)
% Default is [1 2 3] (no permutation)
% 
% OUTPUT)
% C3Ddata : C3Ddata strucure with added joint torques
% 
% NOTES)
% This file might require expansion for other joint torques

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
fid = fopen(subjosidsto,'r');
if fid == -1
    warning(['Failed to open ' subjosidsto '. Skipping.']);
    return;
end

tline = '';
ntimeout = 0;
while isempty(strfind(tline,'time')) && ntimeout < 20
    
    tline = fgetl(fid);
    
    if ~isempty(strfind(tline,'nRows='))
        nRows = str2double( tline(strfind(tline,'nRows=') + 6 : end) );
    elseif ~isempty(strfind(tline,'nColumns='))
        nCols = str2double( tline(strfind(tline,'nColumns=') + 9 : end) );
    end
    
    ntimeout = ntimeout + 1;
end
fclose(fid);

if ntimeout == 20
    error('getOSIDinM:DataHeader','Unable to find data column header in file.');
else
    trqData = dlmread(subjosidsto,'\t',[ntimeout,0,ntimeout+nRows-1,nCols-1]);
    % NOTE: assumed here the data is directly below the column header
end


%% Sort data based on header
% 
% NOTE: It is assumed the column header is delimited with tabs
% NOTE2: Even though the torques are not in global coordinates, they are
% sorted in a "global as possible" way. For example, the knee
% flexion-extension moment is assumed to be in the sagittal plane, about
% the transversal axis, regardless of e.g. the upper leg rotation.
% As a reminder, in OpenSim, x is the walking direction, y is upward and z
% to the right

headercell = strsplit(tline,'\t');

% % Column and layer index for data matrix
% stoidx = zeros(length(headercell),3);
% 
% % Create indices and labelcells for data storage based on header names
% for ihd = 1:length(headercell)
%     switch headercell{ihd}
%         case 'pelvis_tx_force'
%             stoidx(ihd,:) = [1 1 ihd];
%         case 'pelvis_ty_force'
%             stoidx(ihd,:) = [1 2 ihd];
%         case 'pelvis_tz_force'
%             stoidx(ihd,:) = [1 3 ihd];
%         case 'pelvis_list_moment'
%             stoidx(ihd,:) = [2 1 ihd];
%         case 'pelvis_rotation_moment'
%             stoidx(ihd,:) = [2 2 ihd];
%         case 'pelvis_tilt_moment'
%             stoidx(ihd,:) = [2 3 ihd];
%           
%         case 'lumbar_bending_moment'
%             stoidx(ihd,:) = [3 1 ihd];
%         case 'lumbar_rotation_moment'
%             stoidx(ihd,:) = [3 2 ihd];
%         case 'lumbar_extension_moment'
%             stoidx(ihd,:) = [3 3 ihd];
%         
%         case 'hip_adduction_l_moment'
%             stoidx(ihd,:) = [4 1 ihd];
%         case 'hip_rotation_l_moment'
%             stoidx(ihd,:) = [4 2 ihd];
%         case 'hip_flexion_l_moment'
%             stoidx(ihd,:) = [4 3 ihd];
%         case 'knee_angle_l_moment'
%             stoidx(ihd,:) = [5 3 ihd];
%         case 'subtalar_angle_l_moment'
%             stoidx(ihd,:) = [6 1 ihd];
%         case 'ankle_angle_l_moment'
%             stoidx(ihd,:) = [6 3 ihd];
%         
%         case 'hip_adduction_r_moment'
%             stoidx(ihd,:) = [7 1 ihd];
%         case 'hip_rotation_r_moment'
%             stoidx(ihd,:) = [7 2 ihd];
%         case 'hip_flexion_r_moment'
%             stoidx(ihd,:) = [7 3 ihd];
%         case 'knee_angle_r_moment'
%             stoidx(ihd,:) = [8 3 ihd];
%         case 'subtalar_angle_r_moment'
%             stoidx(ihd,:) = [9 1 ihd];
%         case 'ankle_angle_r_moment'
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
% labelcellAll = {'PELVIST','PELVISM','LJC','HJCL','KJCL','AJCL','HJCR','KJCR','AJCR'};
% 
% iCols = unique(stoidx(:,1));
% jointTrqDataLabel = labelcellAll(iCols);
% 
% 
% %% Assign data
% 
% % Make data equally long as the rest of the marker data using the sync indices
% msiz = size(C3Ddata.Resample.IMU,1);
% 
% % Pre-alloc new data container
% jointTrqData = zeros(msiz,numel(iCols),3);
% 
% % Put the data in a new matrix
% for ihd = 1:size(stoidx,1)
%     jointTrqData(:,stoidx(ihd,1),stoidx(ihd,2)) = trqData(:,stoidx(ihd,3));
% end
% 
% % Permutate the data dimensions
% jointTrqData = jointTrqData(:,:,permvec);

% Add data to C3DProc structure
C3Ddata.Resample.Sych.IDTrqData_Portable = trqData;
C3Ddata.Resample.Sych.IDTrqDataLabel_Portable = headercell;

end