function [C3Ddata] = getProbeGlobal(C3Ddata,Body)
%% getProbeGlobal
% Reconstruct the global coordinates of a probed position in a measurement
%
% INPUT)
% C3Ddata: data structure of a probe trial containing at least the fields:
%     Marker.MarkerData
%     Marker.MarkerDataLabel
% Body: data structure with the fields:
%     ProbePos{i}.PosName :     name of each probe position i
%     ProbePos{i}.ClusterName : name of each cluster the probe position belongs to
%     ProbePos{i}.UsedLEDs :    LEDs used within the cluster (e.g. cluster has 4 LEDs, but only 3 were visible during probing)
%     ProbePos{i}.PosInLocal
% 
% OUTPUT)
% C3Ddata: same as input, but with added fields:
%     Marker.ProbedData : probe positions
%     Marker.ProbedDataLabel : probe positions labels
% 
% NOTES)
% Requires soder.m 
% 
% Update history)
% 14-11-2014 : Mark Vlutters

%% Reconstruction

    % Construct the probe position (this will take a while)
%     h_wait = waitbar(0);  % Fancy waitbar
    for iPos = 1:length(Body.ProbePos)
        
        % Get information from Body
        posName = Body.ProbePos{iPos}.PosName;
        clusterLabel = Body.ProbePos{iPos}.ClusterLabel;
        usedLEDs = Body.ProbePos{iPos}.UsedLEDs;
        clusterMarkerLocal = Body.ProbePos{iPos}.ClusterMarkerLocal;
        
        MRKR = find( strncmpi(clusterLabel,C3Ddata.Marker.MarkerDataLabel,length(clusterLabel)) );
        MRKR = MRKR(usedLEDs);
        clusterMarkerGlobal = C3Ddata.Marker.MarkerData(:,MRKR,:);
        
        % Remove invisible markers (only when using 4 or more per cluster, otherwise try reconstruct anyway)
        % And remove markers that show too much relative movement to other markers within the cluster
        interTol = 0.01; % Tolerance for shortest interdistance (in meter) : not taken "too" small
        if size(clusterMarkerGlobal,2) > 3
            
            maskVisible = sum(clusterMarkerGlobal,3) ~= 0;
            
            % This is a bit double, as invisible markers generally don't have appropriate inter-distance either
            nWithinTol = zeros(size(clusterMarkerGlobal,1),size(clusterMarkerGlobal,2));
            for imark = 1:size(clusterMarkerGlobal,2)
                interDist = sqrt(sum( ( clusterMarkerGlobal - repmat( clusterMarkerGlobal(:,imark,:) ,[1 size(clusterMarkerGlobal,2) 1]) ).^2 ,3));
                nWithinTol(:,imark) = sum( abs( interDist - repmat(median(interDist,1) , [size(interDist,1) 1]) ) < interTol , 2)  -1 ; % Number of OTHER LEDs (hence -1) that are within tol
            end

            % Remove a marker completely if it has too little useful contribution (less than 10%)
            % Unless you would end up with less than 3 LEDs
            prct = 0.1;
            sumWithinTol = sum( nWithinTol >=2 ,1);
            nWithinTol(:,sumWithinTol < prct.*size(nWithinTol,1)) = 0;
            
            % Check where 3 or more useful LEDs remain
            hasThree = sum( nWithinTol >2 ,2 ) >= 3;
            
            % If there are not >=3 markers sufficing the constraint, select those that are probably the "best" 3
            sortLEDs = sort( sumWithinTol , 'descend');
            bestThree = sortLEDs(1:3);
            [~,idxIn,~] = intersect(sumWithinTol,bestThree);
            nWithinTol(~hasThree,idxIn) = 2;
            
            % NOTE: this can still go wrong, e.g. if you have 6 LEDs for a single cluster, and they somehow "split" in 2 groups of 3.
            maskTol = nWithinTol >= 2;
            
            usedMrk = (maskVisible + maskTol) == 2;
            
        else 
            % There's nothing you can do about bad data now, you cannot reconstruct with less than 3 LEDs
            usedMrk = true(size(clusterMarkerGlobal,1),size(clusterMarkerGlobal,2));
        end
        
        % Reconstruct
        dumVar = zeros(size(clusterMarkerGlobal,1),3);
        for idx = 1:size(C3Ddata.Marker.MarkerData,1)

            [~,d,~] = soder(clusterMarkerLocal(usedMrk(idx,:),:),squeeze(clusterMarkerGlobal(idx,usedMrk(idx,:),:)));
            dumVar(idx,:) = d;
%             eval([posName '(idx,:) = d ;']); % Beware : dynamic variable creation

%             waitbar(idx/size(C3Ddata.Marker.MarkerData,1),h_wait,['Constructing ' posName]);
        end

        % Add to processed marker data
        if ~isfield(C3Ddata.Marker,'ProbedData')
%             eval(['C3Ddata.Marker.ProbedData = zeros(size(' posName ',1),length(Body.ProbePos),3);']);
            C3Ddata.Marker.ProbedData = zeros(size(dumVar,1),length(Body.ProbePos),3);
        end
%         eval(['C3Ddata.Marker.ProbedData(1:length(' posName '),iPos,1:3) = ' posName ';']);
        C3Ddata.Marker.ProbedData(1:size(dumVar,1),iPos,1:3) = dumVar;

        C3Ddata.Marker.ProbedDataLabel{iPos} = posName;
        C3Ddata.Marker.ProbedDataLabelC{iPos} = clusterLabel; % Allows you to trace back probe positions to their corresponding clusters
        
        C3Ddata.Marker = orderfields(C3Ddata.Marker);
        C3Ddata.Analog = orderfields(C3Ddata.Analog);

        % Clean up the created variable
        eval(['clear ' posName]);

    end
%     close(h_wait);

end