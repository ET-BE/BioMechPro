function gBMP_modSort(hObject,eventdata,h1338)
% Sort module tabs

% Get module handles
htabgMod = findobj(h1338,'tag','tabgMod');
htabs = get(htabgMod,'children');

if ~isempty(htabs)

    % Get all tab names
    tabnames = cell(1,length(htabs));
    for itab = 1:length(htabs)
       tabnames{itab} = get(htabs(itab),'title');
    end
    
    % Ask what to do
    butout = questdlg('Move module tabs how?','Module tab sort','Alphabetically','Move current','Move current');


    if strcmp(butout,'Alphabetically')
        
        % Sort alphabetically
        [~,perm] = sort(tabnames);
        htabgMod.Children = htabgMod.Children(perm);

    elseif strcmp(butout,'Move current')
        
        % Get current tab index in the child list
        ctabtitle = get(htabgMod.SelectedTab,'title');
        curidx = find(strcmpi(ctabtitle,tabnames),1,'first');

        % Ask to which position to switch (BEFORE selected)
        [newidx,ok] = listdlg('liststring',tabnames,...
        'name','ModSelect',...
        'promptstring','Select position to move to:',...
        'selectionmode','single');
        
        % Return if cancelled
        if (~ok) || (newidx == curidx)
            return;
        end
        
        % Swap tabs
        if curidx > newidx
            perm = [1:newidx-1 curidx newidx:curidx-1 curidx+1:length(htabgMod.Children)];
        else
            perm = [1:curidx-1 curidx+1:newidx curidx newidx+1:length(htabgMod.Children)];
        end
        
        htabgMod.Children = htabgMod.Children(perm);
        
    else
        return;
    end    
else
    return;
end


end