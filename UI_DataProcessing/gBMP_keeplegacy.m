function [] = gBMP_keeplegacy(hObject,eventdata)
% Check or uncheck keep legacy files
    
    if strcmpi( get(hObject,'checked') , 'off');
        set(hObject,'checked','on');
    elseif strcmpi( get(hObject,'checked'), 'on');
        set(hObject,'checked','off');
    end

end