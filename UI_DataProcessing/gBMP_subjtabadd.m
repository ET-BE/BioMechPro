function [] = gBMP_subjtabadd(hObject,eventdata,h1338)
% Add a tab to the subject tab

% Get subject tabgroup handle
tabgSubj = findobj(h1338,'tag','tabgSubj');

% Count the number of existing tabs
nsubjTabs = length(get(tabgSubj,'children')) - 2; % -2 for + and - tabs

if nsubjTabs >= 20  % 20 subjects max
    
    set(tabgSubj,'selectedtab',findobj(h1338,'tag',['tabSubj' num2str(nsubjTabs)]));
    beep;
    return;
    
else
    
    % Get handle from last subject tab
    hprevtab = findobj(h1338,'tag',['tabSubj' num2str(nsubjTabs)]);
    
    % Create new tab
    hnewtab = uitab(tabgSubj,...
    'backgroundcolor',get(hprevtab,'backgroundcolor'),...
    'foregroundcolor',get(hprevtab,'foregroundcolor'),...
    'title',['Subj' num2str(nsubjTabs + 1)],...
    'tag',['tabSubj' num2str(nsubjTabs + 1)]);

    % Copy all previous tab properties to the new one
    % WARNING: unless using legacy options, copyobj does not copy callbacks
    copyobj( get(hprevtab,'children') , hnewtab );

    % Edit the copied children
    foo = get(hnewtab,'children');
    for iobj = 1:length(foo)
        
        % Rename all tags
        tagstr = get(foo(iobj),'tag');
        if nsubjTabs < 10
            set(foo(iobj),'tag',[tagstr(1:end-1) num2str(nsubjTabs+1)]);
        else
            set(foo(iobj),'tag',[tagstr(1:end-2) num2str(nsubjTabs+1)]);
        end
        
        % Clear edit fields
        if strcmpi(get(foo(iobj),'style'),'edit')
            set(foo(iobj),'string','');
        end
        
        % Clear checkboxes
        if strcmpi(get(foo(iobj),'style'),'checkbox')
            set(foo(iobj),'value',0);
        end
        
    end
    
    % Set new callback for rootfolder pushbutton
    set( findobj(h1338,'tag',['pshRoot' num2str(nsubjTabs+1)]) , 'callback',{@gBMP_getfile,h1338,['edtRoot' num2str(nsubjTabs+1)]});

    % Set new callback for exclusive checkbox
    set( findobj(h1338,'tag',['chkExclusive' num2str(nsubjTabs+1)]) , 'callback',{@gBMP_subjsetexclusive,h1338} );
    
    % Make new tab the current one
    set(tabgSubj,'selectedtab',hnewtab);
    
end


end