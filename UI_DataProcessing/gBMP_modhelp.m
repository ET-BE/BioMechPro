function [] = gBMP_help(hObject,eventdata,h1338)
% Open module documentation

htab = get(findobj(h1338,'tag','tabgMod'),'selectedtab');

if ~isempty(htab)
    eval(['doc ' get(htab,'title')]);
else
    disp('No module opened: cannot display help');
end

end