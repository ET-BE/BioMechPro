function [clusterLabel,clusterMat] = getClusterLabel(labelcell)
%% checkClusterLabel
% Compute cluster names (if any)
% 
% INPUT)
% labelcell: a labelcell
% 
% OUTPUT)
% clusterLabel: longest possible name shared by 3 or more names in labelcell
% clusterMat: NxM matrix, where N is the number of clusters and M the number of names in each cluster
%         Each element is an index of the label in labelcell.
% 
% NOTE)
% clusterMat can contain zeros, e.g. when some clusters consists of 3 markers and some of 4

%% Do stuff

    clusterMat = [];
    
    n = 0; % Another loop variable
    for i = 1:length(labelcell)

        mask = (i == clusterMat); % Check if you already used the marker in a cluster
        if isempty(mask) || ~any(mask(:))

            % Break down mystr until it shares an expression with 3 or more others
            mystrfull = labelcell{i};
            j = length(mystrfull); nmatch = 0;
            while (j>0) && (sum(nmatch) < 3)
                mystr = mystrfull(1:j);
                nmatch = strncmp(mystr,labelcell,length(mystr));
                j = j - 1;
            end

            if sum(nmatch) >= 3
                n = n + 1;
                clusterLabel{n} = mystr;  % If you used proper names, these should correspond with the cluster names.
                clusterMat(n,1:length(find(nmatch))) = find(nmatch);
            end
            
        end
    end


end