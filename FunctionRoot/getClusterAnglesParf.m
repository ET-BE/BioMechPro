function [C3Ddata] = getClusterAnglesParf(C3Ddata,outMode)
%% getSegmentAngles
% Get segment rotation in Euler angles, with respect to the static
% Rotation angles are in radian, and assumed rotated in order Rx, Ry, Rz
% 
% INPUT)
% C3Ddata with at least the fields:
%     .Marker.MarkerData
%     .Marker.MarkerDataStatic
%     .Marker.MarkerDataLabel (required if segmat is not supplied)
% 
% outMode
% 
% OUTPUT)
% C3Ddata with added field(s):
%     C3Ddata.Marker.SegAngle : MxNx3 matrix (data,segments,dimensions) 
%     C3Ddata.Marker.SegLabel : Labels of the N segments. Will be empty if segmat is provided.
% 
% NOTES)
% Requires multiple processor cores to be available
% 
% This requires all markers on the same segment to share a name (e.g. LFOOT1, LFOOT2, LFOOT3).
% It is also assumed that a marker cannot be part of 2 or more clusters.
% 
% Requires getClusterLabel
% Requires soder.m 
% Requires R2EA.m if you want to output Euler Angles
% 
% It might be convenient to have clusterMat as custom input if you don't want
% to use the labels for detecting clusters.

%% Check input

if ~isempty(outMode)
    if ~any( strcmpi(outMode,{'E','R','ER','RE'}) )
        warning('outMode not recognized, setting to Euler angles');
        outMode = 'E';
    end
    if strcmpi(outMode,'RE')
        outMode = 'ER';
    end
else
    outMode = 'E';
end

% if ~isfield(C3Ddata.Marker,'MarkerDataStatic')
%     warning('No static positions present. Taking arbitrary reference (median first 10 samples).');
% end

%% Do stuff

% get clusterLabel and indices of each cluster in the total marker data
[clusterLabel,clusterMat] = getClusterLabel(C3Ddata.Marker.MarkerDataLabel);

% Preallocate segment angle data
if any( strcmpi(outMode,{'E','ER'}) )
    clusterAngleE = zeros(size(C3Ddata.Marker.MarkerData,1) , size(clusterMat,1) , 3 );
end
if any( strcmpi(outMode,{'R','ER'}) )
    clusterAngleR = zeros(size(C3Ddata.Marker.MarkerData,1) , size(clusterMat,1) , 3 , 3);
end

% Construct cluster angles w.r.t. the static
% h_wait = waitbar(0);  % Fancy waitbar
for i = 1:size(clusterMat,1)
    
    % Marker indices that make up the cluster (because clusterMat can contain zeros)
    clusidx = setdiff( clusterMat(i,:) , 0 );
    
%     waitbar(i/size(clusterMat,1),h_wait,['Constructing angles of cluster ' clusterLabel{i}]);
    
    statData = reshape( C3Ddata.Marker.MarkerDataStatic(:,clusidx,:) , [length(clusidx) 3] );
    markData = C3Ddata.Marker.MarkerData(:,clusidx,:);
    
    if strcmpi(outMode,'E')
        
        dumVar1 = zeros(size(markData,1),3);
        marksiz = size(markData); marksiz2 = marksiz(2); marksiz3 = marksiz(3);
        parfor j = 1:size(C3Ddata.Marker.MarkerData,1)
            % Find rotation angle with respect to static posture
            [R,~,~] = soder( statData , reshape(markData(j,:,:), [marksiz2 marksiz3] ) ); 

            dum = R2EA(R);
            dumVar1(j,:) = dum(:,1);
        end 
        clusterAngleE(:,i,1:3) = dumVar1; % ccw
        
    elseif strcmpi(outMode,'R')
        
        dumVar2 = zeros(size(markData,1),3,3);
        marksiz = size(markData); marksiz2 = marksiz(2); marksiz3 = marksiz(3);
        parfor j = 1:size(C3Ddata.Marker.MarkerData,1)
            % Find rotation angle with respect to static posture
            [R,~,~] = soder( statData , reshape(markData(j,:,:), [marksiz2 marksiz3] ) ); 
            dumVar2(j,:,:) = R;
        end 
        clusterAngleR(:,i,1:3,1:3) = dumVar2; % as rotation matrix
        
    else
        
        dumVar1 = zeros(size(markData,1),3);
        dumVar2 = zeros(size(markData,1),3,3);
        marksiz = size(markData); marksiz2 = marksiz(2); marksiz3 = marksiz(3);
        parfor j = 1:size(C3Ddata.Marker.MarkerData,1)
            % Find rotation angle with respect to static posture
            [R,~,~] = soder( statData , reshape(markData(j,:,:), [marksiz2 marksiz3] ) ); 

            dum = R2EA(R);
            dumVar1(j,:) = dum(:,1);
            dumVar2(j,:,:) = R;
        end 
        clusterAngleE(:,i,1:3) = dumVar1; % ccw
        clusterAngleR(:,i,1:3,1:3) = dumVar2; % as rotation matrix
        
    end
    
    
end
% close(h_wait);

if any( strcmpi(outMode,{'E','ER'}) )
    C3Ddata.Marker.ClusterAngleE = clusterAngleE; % Euler format
end
if any( strcmpi(outMode,{'R','ER'}) )
    C3Ddata.Marker.ClusterAngleR = clusterAngleR; % Rotation matrix format
end
C3Ddata.Marker.ClusterLabel = clusterLabel; % The clusters on which you have rotation info

end