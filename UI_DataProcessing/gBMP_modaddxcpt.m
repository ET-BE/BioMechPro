function [] = gBMP_modaddxcpt(hObject,eventdata,h1338)
% Add subject exception tab corresponding to module

% Selected module tab handle
htabmod = get( findobj(h1338,'tag','tabgMod') , 'selectedtab');

% Determine module tab index
modidx = get(htabmod,'tag');
unders = strfind(modidx,'_');
modidx = modidx(unders(end)+1:end);

% Exception tabgroup handle of active module tab
htabgexc = findobj(h1338,'tag',['tabgExc_' modidx]);

% Check for existing exception tabs for active module
tabExcAll = get(htabgexc,'children');
subjExist = zeros(1,length(tabExcAll));
if ~isempty( tabExcAll )
    for itab = 1:length(tabExcAll)
        % Get subjects numbers for which exception already exists
        tabtag = get( tabExcAll(itab) , 'tag');
        unders = strfind(tabtag,'_');
        subjExist(itab) = str2double(tabtag(unders(end)+1:end));
    end
end

% Get the number of subjects (regardless of ignore or exclusive)
nsubjtab = length(get(findobj(h1338,'tag','tabgSubj'),'children'))-2; % -2 for +- tab correct

% Remove subjects for which exception tab already exists
subjInclude = setxor(1:nsubjtab,subjExist);

% Ask user for which subject(s) an exception tab should be made
listcell = cell(1,numel(subjInclude));
for istr = 1:numel(subjInclude);
    listcell{istr} = num2str(subjInclude(istr));
end
[idxsel,ok] = listdlg('liststring',listcell,...
    'name','SubjSelect',...
    'promptstring','Select subject(s):',...
    'selectionmode','multiple');

% Return if cancelled
if ~ok
    return;
end

% Get all module tab children handles (tMod covers both edtMod and txtMod)
hcopy = findobj(htabmod,'-regexp','tag','tMod');

for idxsel = idxsel

    % Create new exception tabs (note: listcell{idxsel} = isubj)
    hnewtab = uitab(htabgexc,...
    'backgroundcolor',[0.85 0.85 0.9],...
    'foregroundcolor',[0 0 0],...
    'title',['SubjExc' listcell{idxsel}],...
    'tag',['tabExc_' num2str(modidx) '_' listcell{idxsel}]);

    % Copy all to new exception tab
    copyobj(hcopy,hnewtab);
    
    % Adjust the text object background color for pretty purposes
    htxtall = findobj(hnewtab,'-regexp','tag','txt');
    for itxt = 1:length(htxtall)
        set(htxtall(itxt),'backgroundcolor',get(hnewtab,'backgroundcolor'));
    end
    
end

% Switch to last created tab
set(htabgexc,'selectedtab',hnewtab);

% Make tabgroup visisble
set(htabgexc,'visible','on');

% Enable remove exception button
set(findobj(h1338,'tag','psh&-'),'enable','on');

end