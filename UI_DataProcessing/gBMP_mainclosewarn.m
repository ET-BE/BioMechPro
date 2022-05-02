function [] = gBMP_mainclosewarn(hObject,eventdata,h1338)
% Warning dialog when closing the UI

    dlgans = questdlg({'Are you sure you want to quit?','(Unsaved data will be lost).'},...
        'Quit?',...
        'No','Yes','No');

    if strcmp(dlgans,'Yes')
        delete(h1338);
    end
    
end