function [C3Ddata] = getClusterAngles(C3Ddata)
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
% OUTPUT)
% C3Ddata with added field(s):
%     C3Ddata.Marker.SegAngle : MxNx3 matrix (data,segments,dimensions) 
%     C3Ddata.Marker.SegLabel : Labels of the N segments. Will be empty if segmat is provided.
% 
% NOTES)
% This requires all markers on the same segment to share a name (e.g. LFOOT1, LFOOT2, LFOOT3).
% It is also assumed that a marker cannot be part of 2 or more clusters.
% 
% Requires getClusterLabel
% Requires soder.m 
% Requires R2EA.m if you want to output Euler Angles
% 
% It might be convenient to have clusterMat as custom input if you don't want
% to use the labels for detecting clusters.

%% Do stuff

% get clusterLabel and indices of each cluster in the total marker data
[clusterLabel,clusterMat] = getClusterLabel(C3Ddata.Marker.MarkerDataLabel);

% Preallocate segment angle data
clusterAngleE = zeros(size(C3Ddata.Marker.MarkerData,1) , size(clusterMat,1) , 3 );
clusterAngleR = zeros(size(C3Ddata.Marker.MarkerData,1) , size(clusterMat,1) , 3 , 3);

% Construct cluster angles w.r.t. the static
h_wait = waitbar(0);  % Fancy waitbar
for i = 1:size(clusterMat,1)
    
    % Logic array of markers that make up the cluster
    clusidx = setdiff( clusterMat(i,:) , 0 );
    
    for j = 1:size(C3Ddata.Marker.MarkerData,1)
    
        % Find rotation angle with respect to static posture
        [R,~,~] = soder( squeeze( C3Ddata.Marker.MarkerDataStatic(:,clusidx,:) ) , squeeze(C3Ddata.Marker.MarkerData(j,clusidx,:)) ); 

        dum = R2EA(R);
        clusterAngleE(j,i,1:3) = dum(:,1); % ccw
        clusterAngleR(j,i,1:3,1:3) = R; % as rotation matrix
        
        waitbar(j/size(C3Ddata.Marker.MarkerData,1),h_wait,['Constructing angles of cluster ' clusterLabel{i}]);
    end    
end
close(h_wait);

C3Ddata.Marker.ClusterAngleE = clusterAngleE; % Euler format
C3Ddata.Marker.ClusterAngleR = clusterAngleR; % Rotation matrix format
C3Ddata.Marker.ClusterLabel = clusterLabel; % The clusters on which you have rotation info

end