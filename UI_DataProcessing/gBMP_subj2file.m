function [] = gBMP_subj2file(hObject,eventdata,h1338)
% Write subject info from the ui to a .m file

% uiputfile
[filename,pathname] = uiputfile('.m');

fullpath = [pathname filename];

if sum(fullpath) ~= 0

    fid = fopen(fullpath,'w');
    if fid == -1
        error('gBMP_subj2file:noopen','Failed to open file to write');
    end
    
    % Write switch statement
    fprintf(fid,'switch isubj \n');
    
    % Get the number of subjects
    tabgSubj = findobj(h1338,'tag','tabgSubj');
    tabsSubj = get(tabgSubj,'children');
    nsubj = length(tabsSubj) - 2;
    
    for isubj = 1:nsubj
                
        % Get info from the fields from subject tabgroup 
        % (+2 to correct for + and - tabs)
        rootrhs = get( findobj(tabsSubj(isubj+2),'-regexp','tag','edtRoot') , 'string');
        trialrhs = get( findobj(tabsSubj(isubj+2),'-regexp','tag','edtTrials') , 'string');
        otherstr = get( findobj(tabsSubj(isubj+2),'-regexp','tag','edtOther') , 'string');
        
        try
            % Write info to file
            fprintf(fid,['\t case ' num2str(isubj) '\n']);
            if ~isempty(rootrhs)
                fprintf(fid,'\t\t subjroot = ''%s''; \n',rootrhs);
            end
            if ~isempty(trialrhs)
                fprintf(fid,'\t\t subjtrials = %s \n',trialrhs);
            end
            if ~isempty(otherstr)
                for iline = 1:size(otherstr,1)
                    fprintf(fid,'\t\t %s \n',otherstr(iline,:));
                end
            end
        catch
            fclose(fid);
            warning('gBMP_subj2file:noclose','Error in writing to file.');
            return;
        end
        
    end
    
    fprintf(fid,'end');
    fclose(fid);
    
else
    return;
end


end