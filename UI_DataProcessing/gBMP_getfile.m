function [] = gBMP_getfile(hObject,eventdata,h1338,tag2put)
% Get path and filename and put in specified edit field

% Select folder
foldername = uigetdir(userpath,'Select subject root folder');

if foldername
    % Put folder in edit field with given tag
    set(findobj(h1338,'tag',tag2put),'string',foldername);
end

end