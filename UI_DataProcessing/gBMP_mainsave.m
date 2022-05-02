function [] = gBMP_mainsave(hObject,eventdata,h1338)
% Save information in UI editable fields
% User must choose a file name if the file does not exist

% Find every child that is an editable field or checkbox
childdata = allchild(h1338);

% Determine whether to use save or save as
if ~exist(get(h1338,'userdata'),'file')
    
    % uiputfile
    [filename,pathname] = uiputfile('.mat');
    
    fullpath = [pathname filename];
    
else
    
    fullpath = get(h1338,'userdata');
    
end

if sum(fullpath) ~= 0

    % Set last path as active file
    set(h1338,'userdata',fullpath);
    
    % Save (cannot input path with savefig, but save gives annoying warning
    % for saving figure handles: only compatible in 2014b and later).
    warning('off','all');
    save(fullpath,'childdata');
    warning('on','all');

else
    return;
end

end