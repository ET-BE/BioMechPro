function [] = gBMP_mainopen(hObject,eventdata,h1338)
% Load previously saved UI subject and module info

    % Ask user to select file
    [filename,pathname] = uigetfile('.mat');

    fullpath = [pathname filename];
    
    
    if sum(fullpath) ~= 0
        % Destroy the current UI
        delete(h1338);
    
        % Load new UI
        load(fullpath);
    end

end