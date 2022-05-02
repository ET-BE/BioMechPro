function [] = gBMP_modcloseall(hObject,eventdata,h1338)
% Close all modules

% Find module tabgroup handle
htabgmod = findobj(h1338,'tag','tabgMod');

% Return if no modules exist
if isempty(get(htabgmod,'children'))
    return;
end

% Ask for closure
dlgans = questdlg({'Are you sure you want to close all modules?','(Unsaved data will be lost).'},...
        'Close all?',...
        'No','Yes','No');

if ~strcmp(dlgans,'Yes')
    return;
end

% Find all exception tabgroup handles
htabgExcAll = findobj(h1338,'-regexp','tag','tabgExc');

% Destroy all children (-insert evil laugh here-)
delete( htabgExcAll );
delete( get(htabgmod,'children') );

% Make module tabgroup invisible and disable the module buttons
if isempty( get(findobj(h1338,'tag','tabgMod'),'selectedtab') )
    set(htabgmod,'visible','off');
    set(findobj(h1338,'tag','psh>'),'enable','off');
    set(findobj(h1338,'tag','psh&+'),'enable','off');
    set(findobj(h1338,'tag','psh&-'),'enable','off');
end


end