function [copdata] = getCOP(fdata,fthresh)
% Get COP data
% 
% INPUT)
% fdata : Nx6 single plate or Nx12 dual plate data
% fthresh : force threshold in N for determining contact
% 
% OUTPUT)
% copdata : structure with 1 or 3 Nx3 matrices, each containing COPx, COPy, COPz
% In case of 3 matrices, it contains the COP for left, right and combined plates
% Data below fthresh is made NaN.
% x -> AP
% y -> vert
% z -> ML
%
% NOTES)
% Input is assumed:
% an Nx12 matrix with
% 1-3 : Fx, Fy, Fz
% 4-6 : Mx, My, Mz
% 
% OR
% 
% an Nx12 matrix with : 
% 1-3 : Fx, Fy, Fz of the left plate
% 4-6 : Mx, My, Mz of the left plate
% 7-9 : Fx, Fy, Fz of the right plate
% 10-12 : Mx, My, Mz of the right plate
% 
% In a right-handed coordinate system, where z is the vertical axis
% The function assumes that COPz (vertical) is zero
% 
% The total COP of two plates only has a meaning if both plates have the same origin.

%% Input data check
    if size(fdata,2) == 6
        nplate = 1;
    elseif size(fdata,2) == 12
        nplate = 2;
    else
        error('getCOP:input data does not have 6 or 12 columns');
    end


%% Get COP
% COPz is assumed zero

    if nplate == 1
        % Set contactless force to NaN
        fdata( fdata(:,2) < fthresh , 2) = NaN;
        
        % Get COP
        COPT = [fdata(:,6) ./ fdata(:,2) , zeros(length(fdata),1), -fdata(:,4) ./ fdata(:,2)];

        % Store
        copdata.COPT = COPT;
        
    else
        % Fuse both plates for total COP
        fdatat = fdata(:,1:6) + fdata(:,7:12);
        
        % Set contactless force to NaN
        fdata( fdata(:,2) < fthresh , 2) = NaN;
        fdata( fdata(:,8) < fthresh , 8) = NaN;
        
        % Get COP
        COPL = [ fdata(:,6) ./ fdata(:,2) , zeros(length(fdata),1), -fdata(:,4) ./ fdata(:,2)];
        COPR = [ fdata(:,12) ./ fdata(:,8) , zeros(length(fdata),1), -fdata(:,10) ./ fdata(:,8) ];
        COPT = [ fdatat(:,6) ./ fdatat(:,2) , zeros(length(fdatat),1), -fdatat(:,4) ./ fdatat(:,2)];

        % Store
        copdata.COPL = COPL;
        copdata.COPR = COPR;
        copdata.COPT = COPT;
    end

end
