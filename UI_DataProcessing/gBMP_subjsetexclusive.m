function [] = gBMP_subjsetexclusive(hObject,eventdata,h1338)
% Checkbox to make a subject exclusive for processing, ignoring all others

% Get checkbox info
chkval = get(hObject,'value');
chktag = get(hObject,'tag');

% Current subject tab 
% (This is a dirty workaround to find out if the current subject is <10 
% or not. It stops to function if you change the tag of the checkbox.)
if length(get(hObject,'tag')) == 13
    subjnum = chktag(end);
elseif length(get(hObject,'tag')) == 14
    subjnum = chktag(end-1:end);
end

% Get number of subject tabs
nsubjTabs = length(get(findobj(h1338,'tag','tabgSubj'),'children')) - 2; % -2 for + and - tabs

if chkval
    
    % Disable ignore check of the same subject
%     set( findobj(h1338,'tag',['chkIgnore' chktag(end)]) , 'value',0);
    
    % Disable all other exclusive checks
    for isubj = setdiff( 1:nsubjTabs , str2double(subjnum) )
        set( findobj(h1338,'tag',['chkExclusive' num2str(isubj)]), 'value',0 );
    end
    
end

end