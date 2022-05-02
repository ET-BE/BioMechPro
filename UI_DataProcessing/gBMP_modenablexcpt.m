function [] = gBMP_modenablexcpt(hObject,eventdata,h1338)
% Enable / disable 'exception' buttons and tabgroup depending on selected module tab

% Set all exception tabs visibility to off
htabgexcAll = findobj(h1338,'-regexp','tag','tabgExc');
for itabg = 1:length(htabgexcAll)
    set(htabgexcAll(itabg),'visible','off');
end

% Selected module tab handle
htabmod = get( findobj(h1338,'tag','tabgMod') , 'selectedtab');

% Selected module index
modidx = get(htabmod,'tag');
unders = strfind(modidx,'_');
modidx = modidx(unders(end)+1:end);

% Check if there are any edit fields, and enable things if so
if isempty( findobj(htabmod,'-regexp','tag','edtMod') )
    set(findobj(h1338,'tag','psh&+'),'enable','off');
else
    set(findobj(h1338,'tag','psh&+'),'enable','on');
end

% Check if any exceptions exist within the module.
% If so, make tabgroup visisble and enable remove button.
if ~isempty( get(findobj(h1338,'tag',['tabgExc_' modidx]),'children') );
    set(findobj(h1338,'tag',['tabgExc_' modidx]),'visible','on');
    
    set(findobj(h1338,'tag','psh&-'),'enable','on');
else
    set(findobj(h1338,'tag','psh&-'),'enable','off');
end

end