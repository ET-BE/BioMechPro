function [C3Ddata] = getStatic(C3Ddata,C3Dstatic)
%% getStatic
% Get reference positions from static trial (if any)
% 
% INPUT)
% C3Ddata: the general data structure, with at least the field:
%     .Marker.MarkerDataLabel
% C3Dstatic: general data structure of the static, with at least the field:
%     .Marker.MarkerData
%     .Marker.MarkerDataLabel
% 
% OUTPUT)
% C3Ddata: thhe general data structure with added field:
%     .Marker.MarkerDataStatic
%     .Marker.ProbedDataStatic
% No corresponding label cells are generated, as these are the same as the
% original label cells for the MarkerData and ProbedData
% 
% NOTES)
% It is assumed that the static data in C3Dstatic contains AT LEAST all 
% markers that are also in C3Ddata.
% 
% Obviously, it is assumed the subject does not move during the static trial
% 
% The static is not squeezed, even though it consists of only a single video frame.
% This way the marker number stays in the second dimension, as it is everywhere.

%% Do stuff

% Remove from C3Dstatic any marker that is not in C3Ddata
statLabelM = C3Dstatic.Marker.MarkerDataLabel;
datLabelM = C3Ddata.Marker.MarkerDataLabel;

idxStatM = zeros(1,length(datLabelM));
for i = 1:length(datLabelM)
    idxStatM(i) = find( strcmpi( datLabelM{i} , statLabelM ) );
end

% Get static (median) and store data
C3Ddata.Marker.MarkerDataStatic = median( C3Dstatic.Marker.MarkerData(:,idxStatM,:) , 1);



% Remove from C3Dstatic any probe position that is not in C3Ddata (if any)
if isfield(C3Dstatic.Marker,'ProbedDataLabel') && isfield(C3Ddata.Marker,'ProbedDataLabel')
    statLabelP = C3Dstatic.Marker.ProbedDataLabel;
    datLabelP = C3Ddata.Marker.ProbedDataLabel;
    
    idxStatP = zeros(1,length(datLabelP));
    for i = 1:length(datLabelP)
        idxStatP(i) = find( strcmpi( datLabelP{i} , statLabelP ) );
    end
    
    % Get static (median) and store data
    C3Ddata.Marker.ProbedDataStatic = median( C3Dstatic.Marker.ProbedData(:,idxStatP,:) , 1);
end

% Order data
C3Ddata.Marker = orderfields(C3Ddata.Marker);

end