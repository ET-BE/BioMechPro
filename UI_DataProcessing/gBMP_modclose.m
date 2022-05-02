function [] = gBMP_modclose(hObject,eventdata,h1338)
% Close module

% Find module tabgroup handle
htabgmod = findobj(h1338,'tag','tabgMod');

% Return if no modules exist
if isempty(get(htabgmod,'children'))
    return;
end

% Ask for closure
dlgans = questdlg({'Are you sure you want to close the module?','(Unsaved data will be lost).'},...
        'Close?',...
        'No','Yes','No');

if ~strcmp(dlgans,'Yes')
    return;
end

% Get handle to selected module and its index
htabmod = get(htabgmod,'selectedtab');
modidx = get(htabmod,'tag');
unders = strfind(modidx,'_');
modidx = modidx(unders(end)+1:end);

% Delete exception tabgroup with same index
delete( findobj(h1338,'tag',['tabgExc_' num2str(modidx)]) );

% Delete active module
delete( htabmod );

% Check if any modules remain, if not, make module tabgroup invisible
% and disable the module buttons
if isempty( get(findobj(h1338,'tag','tabgMod'),'selectedtab') )
    set(htabgmod,'visible','off');
    set(findobj(h1338,'tag','psh>'),'enable','off');
    set(findobj(h1338,'tag','psh&-'),'enable','off');
    set(findobj(h1338,'tag','psh&+'),'enable','off');
else
    % Run the module 'tabchanged' function
    gBMP_modenablexcpt(hObject,eventdata,h1338);
end


end