function [C3Ddata] = correctIKAngle(C3Ddata)
%% correctIKAngle
% Function to correct angles calculated using OpenSim Inverse kinematics
% Sometimes due to "floating" markers or probe points, the IK can assume
% the joint makes a pi or 2pi rad (180 or 360 deg) rotation, which shifts 
% part of the data far away from zero. This function subtracts pi and 2*pi
% from the shifted part(s).
% 
% In general, don't use this function...
 

%% Get data

jointAngData = C3Ddata.Marker.JointAngData;

% Removing 2pi offets
corrbound = 1.2*pi;
corrMask = sum(abs(jointAngData) > corrbound,3) > 0;
corrMaskd = diff(corrMask,1,1);

for icol = 1:size(corrMaskd,2)
    
    corrStart = find(corrMaskd(:,icol) == 1) + 1;
    corrEnd = find(corrMaskd(:,icol) == -1);
    
    % Correct for startless ends and endless starts
    skip = false;
    if ~isempty(corrStart) && ~isempty(corrEnd)
        if corrStart(1) > corrEnd(1)
            corrStart = [1 corrStart];
        end
        if corrStart(end) > corrEnd(end)
            corrEnd = [corrEnd size(corrMask,1)];
        end
    elseif isempty(corrStart) && ~isempty(corrEnd)
        corrStart = 1;
    elseif ~isempty(corrStart) && isempty(corrEnd)
        corrEnd = size(corrMask,1);
    else 
        skip = true;
    end
    
    if ~skip
        for idim = 1:3
            for idx = 1:length(corrStart)
                if (corrEnd(idx) - corrStart(idx)) > 5 % samples
                    if mean(jointAngData(corrStart(idx):corrEnd(idx),icol,idim),1) > corrbound
                        jointAngData(corrStart(idx):corrEnd(idx),icol,idim) = jointAngData(corrStart(idx):corrEnd(idx),icol,idim) - 2*pi;
                    elseif mean(jointAngData(corrStart(idx):corrEnd(idx),icol,idim),1) < -corrbound
                        jointAngData(corrStart(idx):corrEnd(idx),icol,idim) = jointAngData(corrStart(idx):corrEnd(idx),icol,idim) + 2*pi;
                    end
%                     jointAngData(corrStart(idx):corrEnd(idx),icol,idim) = jointAngData(corrStart(idx):corrEnd(idx),icol,idim) - repmat(mean(jointAngData(corrStart(idx):corrEnd(idx),icol,idim),1),[corrEnd(idx)-corrStart(idx)+1 1 1]);
                end
            end
        end
    end
    
end


%% Removing pi offets
corrbound = 0.75*pi;
corrMask = sum(abs(jointAngData) > corrbound ,3) > 0;
corrMaskd = diff(corrMask,1,1);

for icol = 1:size(corrMaskd,2)
    
    corrStart = find(corrMaskd(:,icol) == 1) + 1;
    corrEnd = find(corrMaskd(:,icol) == -1);
    
    % Correct for startless ends and endless starts
    skip = false;
    if ~isempty(corrStart) && ~isempty(corrEnd)
        if corrStart(1) > corrEnd(1)
            corrStart = [1 corrStart];
        end
        if corrStart(end) > corrEnd(end)
            corrEnd = [corrEnd size(corrMask,1)];
        end
    elseif isempty(corrStart) && ~isempty(corrEnd)
        corrStart = 1;
    elseif ~isempty(corrStart) && isempty(corrEnd)
        corrEnd = size(corrMask,1);
    else 
        skip = true;
    end
    
    if ~skip
        for idim = 1:3
            for idx = 1:length(corrStart)
                if (corrEnd(idx) - corrStart(idx)) > 5
                    if mean(jointAngData(corrStart(idx):corrEnd(idx),icol,idim),1) > corrbound
                        jointAngData(corrStart(idx):corrEnd(idx),icol,idim) = jointAngData(corrStart(idx):corrEnd(idx),icol,idim) - pi;
                    elseif mean(jointAngData(corrStart(idx):corrEnd(idx),icol,idim),1) < -corrbound
                        jointAngData(corrStart(idx):corrEnd(idx),icol,idim) = jointAngData(corrStart(idx):corrEnd(idx),icol,idim) + pi;
                    end
%                     jointAngData(corrStart(idx):corrEnd(idx),icol,idim) = jointAngData(corrStart(idx):corrEnd(idx),icol,idim) - repmat(mean(jointAngData(corrStart(idx):corrEnd(idx),icol,idim),1),[corrEnd(idx)-corrStart(idx)+1 1 1]);
                end
            end
        end
    end
    
end


%% Assign data

C3Ddata.Marker.JointAngData = jointAngData;

end