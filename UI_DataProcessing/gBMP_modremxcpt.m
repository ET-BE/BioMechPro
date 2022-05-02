function [] = gBMP_modremxcpt(hObject,eventdata,h1338)
% Remove subject exception from active module
% You don't have to check if exceptions exist, because this function cannot
% be called if no exceptions exist.

% Handle of current module tab
hmodtab = get( findobj(h1338,'tag','tabgMod') , 'selectedtab');

% Get module index
modidx = get(hmodtab,'tag');
unders = strfind(modidx,'_');
modidx = modidx(unders(end)+1:end);

% Handle to exception tabgroup
htabgExc = findobj(h1338,'tag',['tabgExc_' modidx]);

% Find out for which subjects an exception exists
hexcTabAll = get(htabgExc,'children');

subjCell = cell(1,length(hexcTabAll));
for itab = 1:length(hexcTabAll)
    tabtag = get( hexcTabAll(itab) , 'tag');
    unders = strfind(tabtag,'_');
    subjCell{itab} = tabtag(unders(end)+1:end);
end

% Ask user which subjects to remove
[idxsel,ok] = listdlg('liststring',subjCell,...
    'name','SubjSelect',...
    'promptstring','Select subject(s):',...
    'selectionmode','multiple');

if ~ok
    return;
end

% Remove exceptions
for idx = idxsel
    delete( findobj(h1338,'tag',['tabExc_' modidx '_' subjCell{idx}]) );
end

% Check if remain, if not make exception tabgroup invisible and disable
% remove button
if length(idxsel) == length(subjCell)
    set( findobj(h1338,'tag',['tabgExc_' modidx]) , 'visible','off');
    set(findobj(h1338,'tag','psh&-'),'enable','off');
end


end