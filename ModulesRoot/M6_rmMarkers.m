function [Datastr] = M6_rmMarkers(Datastr,field2rm)
% gBMPDynUI field2rm=1;
%
% INPUT)
% - Datastr: the data structure with any number of fields
% 
% - field2rm: string or cell of strings, specifying the field(s) in Datastr
% that should be removed
% 
% - Datastr: the data structure, without the fields specified by field2rm


%% Check 

if ischar(field2rm)
    field2rm = {field2rm};
end

%% Find the fields to remove

% Find
for fld = 1:numel(field2rm)
    
    labelField = field2rm{fld};
    
    layer1 = Datastr.Marker.MarkerDataLabel;
    idx = find(strcmp(layer1, labelField)); % Mind that strcmp is used here (case sensitive)
    if idx==0
        disp([labelField ' not found. Skipping'])
    else
    Datastr.Marker.MarkerDataLabel(idx)=[];
    Datastr.Marker.MarkerData(:,idx,:)=[];
    end
%     disp([labelField ' removed'])
end



end