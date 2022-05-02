function [] = gBMP_mainsaveas(hObject,eventdata,h1338)
% Save information in UI editable fields

% Find every child that is an editable field or checkbox
% TO DO : should pass module type or name for module fields
childdata = allchild(h1338);

% uiputfile
[filename,pathname] = uiputfile('.mat');

fullpath = [pathname filename];

if sum(fullpath) ~= 0

    % Set last path as active file
    set(h1338,'userdata',fullpath);
    
    % Save
    warning('off','all');
    save(fullpath,'childdata');
    warning('on','all');

else
    return;
end

end