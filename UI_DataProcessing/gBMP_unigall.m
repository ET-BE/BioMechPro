function [] = gBMP_unigall(hObject,eventdata,h1338)

% Get subject tabs and number of subjects
tabgSubj = findobj(h1338,'tag','tabgSubj'); 
nsubj = length(get(tabgSubj,'children')) - 2;

% Set all to unignore, unless all are already unignored
igall = true;
for isubj = 1:nsubj
    hchk = findobj(tabgSubj,'tag',['chkIgnore' num2str(isubj)]);
    if get(hchk ,'value')
        set(hchk,'value',0);
        igall = false;
    end
end

if igall
    for isubj = 1:nsubj
        set( findobj(tabgSubj,'tag',['chkIgnore' num2str(isubj)]) , 'value',1);
    end
end


end