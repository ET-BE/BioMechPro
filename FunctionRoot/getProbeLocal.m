function [clusterMarkerLocal,idxC] = getProbeLocal(C3Ddata,cName,idxProbe,probeMarkerLocal)
%% getProbeLocal
% INPUT)
% C3Ddata: data structure of a probe trial containing at least the fields:
%     Marker.MarkerData
%     Marker.MarkerDataLabel
% cName: string, cluster name of the markers that make up the segment corresponding to the
%     probe position file being processed. This name must be present in Marker.MarkerDataLabel.
% idxProbe: vector, indices (columns) of the probe LEDs in Marker.MarkerData.
%     These MUST have the same order as specified in ProbeMarkerLocal.
% ProbeMarkerLocal: matrix, local distances of the probe LEDs to its tip. 
% 
% OUTPUT)
% clusterMarkerLocal: matrix, cluster LED positions with respect to probe position in [0,0,0]
% idxC: vector, indices of cluster LEDs that that provided useful information during probing.
%     These indices are used to select LEDs. For example, you might have a cluster made up 
%     of 4 LEDs, but if only 3 were visible during probing the 4th should 
%     not be used in reconstructing global probe positions later on.
%     Ideally this should have the same value as the number of LEDs your cluster consists of.
% 
% NOTES)
% Requires soder.m !
% 
% Update history)
% 14-11-2014 : Mark Vlutters

%% Find frame distances to probe tip

    % Some numbers that add up to unique numbers
    codevector = [1 2 4 8 16 32];

    % Extract markers of the probe
    probeMarkerGlobal = C3Ddata.Marker.MarkerData(:,idxProbe,:);

    % Data size check
    if size(probeMarkerGlobal,2) ~= size(probeMarkerLocal,1)
        error('getProbeLocal:probeMarkerGlobal does not have the same amount of LEDs as specified in probeMarkerLocal.');
    end
    
    % Extract markers of the relevant cluster (can be more than 3!)
    clusterMarkerGlobal = C3Ddata.Marker.MarkerData(:, strncmpi(cName,C3Ddata.Marker.MarkerDataLabel,length(cName))==1 ,:);
    if size(clusterMarkerGlobal,2) < 3
        error(['getJointsLocal: Insufficient markers with name ' cName],'insufficient markers');
    end
    
    % Data sizes
    datasizP = size(probeMarkerGlobal);
    datasizF = size(clusterMarkerGlobal);
    
    % Check on probe data : isn't a probe led very far away from the rest? (can happen with buggy VZ)
    interTol = 0.01; % Tolerance for shortest interdistance (in meter)
    isWithinTolP = zeros(datasizP(1),datasizP(2));
    for iled = 1:size(probeMarkerLocal,1);
        interLEDGlobalP = sqrt( sum( ( probeMarkerGlobal - repmat(probeMarkerGlobal(:,iled,:),[1 datasizP(2) 1]) ).^2  ,3 ) ); 
        interLEDLocalP = sqrt( sum( ( probeMarkerLocal - repmat(probeMarkerLocal(iled,:),[size(probeMarkerLocal,1) 1]) ).^2 ,2 ) ); 
        isWithinTolP(:,iled) = sum( abs( interLEDGlobalP - repmat(interLEDLocalP' , [size(interLEDGlobalP,1) 1]) ) < interTol ,2) - 1;  % Per LED (col) the number of OTHER LEDs that are within distance of probeMarkerLocal
    end
    
    % Check on cluster data : isn't a cluster LED very far away from the rest?
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
        usedMrk = (maskVisible + maskTol) == 2;  % [nData x nMarker]

        clusterMarkerGlobal( repmat(~usedMrk , [1 1 3]) ) = 0;
            
%     else 
%         % There's nothing you can do about bad data now, you cannot reconstruct with less than 3 LEDs
    end
    
    
    % Find samples where LEDs have less than 2 others nearby, and make these invisible (0)
    probeMarkerGlobal( repmat(isWithinTolP < 2 , [1 1 3]) ) = 0;
    
    % Find as many as possible video frames where THE SAME LEDs are visible, with as many LEDs as possible
    % Note : invisible data = 0
    codemaskP = repmat(codevector(1:datasizP(2)),[datasizP(1),1]); 
    viscodeP = sum( ( sum(probeMarkerGlobal==0,3)==0 ) .* codemaskP ,2);  % LED code per sample
    [~,codesumP] = max( hist(viscodeP,1:sum(codevector)) );               % most prevalent code
    
    codemaskC = repmat(codevector(1:datasizF(2)),[datasizF(1),1]);
    viscodeC = sum( ( sum(clusterMarkerGlobal==0,3)==0 ) .* codemaskC ,2);
    [~,codesumC] = max( hist(viscodeC,1:sum(codevector)) );

    % Find video frames where these LEDs of BOTH PROBE AND FRAME are visible and remove the rest
    idx_vis = ((viscodeP == codesumP) + (viscodeC == codesumC)) == 2;
    probeMarkerGlobal = probeMarkerGlobal( idx_vis , : , :);
    clusterMarkerGlobal = clusterMarkerGlobal( idx_vis , : , :);

    % Decode for markers
    for c = length(codevector):-1:1
        if (codesumP - codevector(c)) >= 0
            codesumP = codesumP - codevector(c);
            codeP(c) = codevector(c);
        else
            codeP(c) = 0;
        end
        if (codesumC - codevector(c)) >= 0
            codesumC = codesumC - codevector(c);
            codeC(c) = codevector(c);
        end
    end

    % Remove markers that are not or little visible throughout the data
    [~,idxP] = intersect(codevector,codeP);
    [~,idxC] = intersect(codevector,codeC);     % You might need this idx in the reconstruction (if >3 frame markers of which 1 or more are unused)!
    probeMarkerGlobal = probeMarkerGlobal(:, idxP , :);
    ProbeMarkerLocal_select = probeMarkerLocal(idxP , :);
    clusterMarkerGlobal = clusterMarkerGlobal(:, idxC , :);

    if size(probeMarkerGlobal ,2) < 3
        error(['getJointsLocal: Probe is poorly or not visible'],'no probe');
    end
    if size(clusterMarkerGlobal ,2) < 3
        error(['getJointsLocal: Marker frame ' cName ' is poorly or not visible'],'no probe');
    end

    % Find 10 indices with as many as possible probe markers visible
%     idx_vis10 = floor( linspace(length(idx_vis)/10,length(idx_vis)-length(idx_vis)/10,10) );  % 10 evenly spaced idx for median value
    idx_vis10 = floor( linspace(size(probeMarkerGlobal,1)/10,size(probeMarkerGlobal,1)-size(probeMarkerGlobal,1)/10,10) );  % 10 evenly spaced idx for median value
    
    clusterMarkerLocal = zeros(size(clusterMarkerGlobal,2),3,length(idx_vis10));
    for igettip = 1:length(idx_vis10)     
        % Check for instance with most markers visible
        idx = idx_vis10(igettip);

        % Determine rotation matrix local probe tip to global probe tip
        [~,d,~] = soder(ProbeMarkerLocal_select,squeeze(probeMarkerGlobal(idx,:,:)));

        % Determine global tip position
        TipGlobal = d'; % is the same as below
    %         TipGlobal = mean( squeeze(ProbeMarkerGlobal(idx,:,:)) - ProbeMarkerLocal*R' , 1);

        % Determine local tip position with respect to relevant frame
        % Contains the probed position w.r.t. the markers in local coordinates
        clusterMarkerLocal(:,:,igettip) = squeeze(clusterMarkerGlobal(idx,:,:)) - repmat(TipGlobal,[size(clusterMarkerGlobal,2),1]);  % nmark,dim,igettip
    end

clusterMarkerLocal = median(clusterMarkerLocal,3);

end