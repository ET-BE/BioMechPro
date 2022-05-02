function [] = gBMP_subjtabrem(hObject,eventdata,h1338)
% Remove a tab from the subject tab

% Get subject tabgroup handle
tabgSubj = findobj(h1338,'tag','tabgSubj');

% Count the number of existing tabs
nsubjTabs = length(get(tabgSubj,'children')) - 2; % -2 for + and - tabs

if nsubjTabs == 1
    
    % Move to last subject tab and do nothing
    set(tabgSubj,'selectedtab',findobj(h1338,'tag',['tabSubj' num2str(nsubjTabs)]));
    beep;
    return;
    
else
    
    % Find the last subject tab and delete it
    delete(findobj(h1338,'tag',['tabSubj' num2str(nsubjTabs)]));
    
    set(tabgSubj,'selectedtab',findobj(h1338,'tag',['tabSubj' num2str(nsubjTabs-1)]));
    
end


end