function [] = gBMP_modrun(hObject,eventdata,h1338)
% Run active module
% Find the module number in UI
% pshtag = get(hObject,'tag');
% modnr = str2double(pshtag(5:end)); % Run tag is psh># , so skip first 4

% Check legacy files toggle and ask for confirmation
if strcmpi( get(findobj(h1338,'tag','uikeeplegacy'),'checked') , 'on')
    questans = questdlg({'Run module?','(WARNING: Legacy option is on.)'},'Confirm','No','Yes','No');
else
    questans = questdlg('Run module?','Confirm','No','Yes','No');
end
if ~strcmp(questans,'Yes')
    return;
end


% Delete the error log file if it exists
% if exist('gBMPErrLog.txt','file')
%     delete('gBMPErrLog.txt');
% end

% Check if legacy files should be maintained
keepLegacy = strcmpi( get(findobj(h1338,'tag','uikeeplegacy') , 'checked') , 'on');

% Find the number of subjects
tabgSubj = findobj(h1338,'tag','tabgSubj');
nsubjTabs = length(get(tabgSubj,'children')) - 2; % -2 for + and - tabs

% Check if a subject is exclusive or ignored
subjects = 1:nsubjTabs;
for isubj = 1:nsubjTabs

    if get( findobj(tabgSubj,'tag',['chkExclusive' num2str(isubj)]) , 'value') % Exclusive overrides Ignore !
        subjects = isubj;
        break;
    elseif get( findobj(tabgSubj,'tag',['chkIgnore' num2str(isubj)]) , 'value')
        subjects(isubj) = NaN;
    end
    
end
subjects(isnan(subjects)) = [];

% Run a check on invalid subject field information (only folder and trials)
for isubj = subjects

    % Check root folder (may only be string)
    subjrootfolder = get(findobj(tabgSubj,'tag',['edtRoot' num2str(isubj)]),'string');
    if isempty(subjrootfolder)
        error('gBMP_modrun:noroot',['Root folder of subject' num2str(isubj) ' is empty']);
    elseif exist(subjrootfolder,'dir') ~= 7
        error('gBMP_modrun:invalidroot',['Invalid root folder for subject ' num2str(isubj)]);
    end
        
    % Check trials (when empty, is considered as ignore)
    subjtrials = regexprep( get(findobj(tabgSubj,'tag',['edtTrials' num2str(isubj)]),'string') ,' ','');
    if isempty(subjtrials)
        warning('gBMP_modrun:emptytrial',['No trials specified for subject ' num2str(isubj) '. Subject will be ignored.']);
%         set(findobj(tabgSubj,'tag',['chkIgnore' num2str(isubj)]) ,'value',1);
        subjects(subjects == isubj) = [];
    else
        try
            otherStr = eval(subjtrials);
            if ( ~((size(otherStr,1)==1)||(size(otherStr,2)==1)) )&&( ~(iscell(subjtrials)||isnumeric(subjtrials)) )
                error('gBMP_modrun:trialvec',['Trials of subject ' num2str(isubj) ' must be a row or column vector or cell']);
            end
        catch
            error('gBMP_modrun:invalidtrials',['Unable to evaluate Trials field for subject ' num2str(isubj)]);
        end
        
    end
    
end

if isempty(subjects)
    warning('gBMP_modrun:nosubj','No subjects to evaluate');
    return;
end

% Determine module tab index
htabmod = get( findobj(h1338,'tag','tabgMod') , 'selectedtab');
modidx = get(htabmod,'tag');
unders = strfind(modidx,'_');
modidx = modidx(unders(end)+1:end);

% Get handle to corresponding exception tabgroup, get exception subject#
htabgexc = findobj(h1338,'tag',['tabgExc_' modidx]);
tabExcAll = get(htabgexc,'children');
subjExc = zeros(1,length(tabExcAll));
if ~isempty( tabExcAll )
    for itab = 1:length(tabExcAll)
        % Get subjects numbers for which exception exists
        tabtag = get( tabExcAll(itab) , 'tag');
        unders = strfind(tabtag,'_');
        subjExc(itab) = str2double(tabtag(unders(end)+1:end));
    end
end

% Check if the module and exception editable fields can be evaluated
hfields = findobj(htabmod,'-regexp','tag','edt');
for ifield = 1:length(hfields)
    fieldstr = get(hfields(ifield),'string');
    
    % Convert to single string (get rid of multiple rows)
    if ~any( size(fieldstr) == 1)
        fieldstr = regexprep(strjoin(cellstr(fieldstr)),'\.\.\.','');
    end
    
    if ~isempty(fieldstr)
        try
            eval(['foo = ' fieldstr ';']);
        catch
            error('gBMP_modrun:eval',['Unable to evaluate module field with content: ' fieldstr]);
        end
    end
end
if ~isempty(subjExc)
    htabexcAll = findobj(htabgexc,'-regexp','tag',['tabExc_' modidx '_']);
    
    for itab = 1:length(htabexcAll)
        htabexc = htabexcAll(itab);
        hfields = findobj(htabexc,'-regexp','tag','edt');
        for ifield = 1:length(hfields)
            fieldstr = get(hfields(ifield),'string');
            
            % Convert to single string (get rid of multiple rows)
            if ~any( size(fieldstr) == 1)
                fieldstr = regexprep(strjoin(cellstr(fieldstr)),'\.\.\.','');
            end
            
            if ~isempty(fieldstr)
                try
                    eval(['foo = ' fieldstr ';']);
                catch
                    error('gBMP_modrun:eval',['Unable to evaluate exception field with content: ' fieldstr]);
                end
            end
        end
    end
    
end

%% Make all edit fields invisible
% Prevent editing / interrupting

set(findobj(h1338,'tag','tabgSubj'),'visible','off');
set(findobj(h1338,'tag','tabgMod'),'visible','off');
set(findobj(h1338,'-regexp','tag','tabgExc'),'visible','off');
set(findobj(h1338,'tag','psh&+'),'enable','off');
set(findobj(h1338,'tag','psh&-'),'enable','off');
set(findobj(h1338,'tag','psh>'),'enable','off');

set(findobj(h1338,'tag','uifile'),'enable','off');
set(findobj(h1338,'tag','uiproc'),'enable','off');
set(findobj(h1338,'tag','uihelp'),'enable','off');

% Give the UI some time to update
pause(0.25)

%% Looping and running the module

% Subject loop
for isubj = subjects

    % Get subject information
    subjrootfolder = get(findobj(tabgSubj,'tag',['edtRoot' num2str(isubj)]),'string');
    subjtrials = eval( get(findobj(tabgSubj,'tag',['edtTrials' num2str(isubj)]),'string') );  % Can be cell or vector!
    
    % If numeric, convert numeric range to cell with strings
    otherStr = cell(1,length(subjtrials));
    if isnumeric(subjtrials)
        for itrial = 1:length(subjtrials)
            otherStr{itrial} = num2str(subjtrials(itrial));
        end
        subjtrials = otherStr;
        clear foo;
    end
    
    % Get all .mat filenames in rootfolder
    rootfiles = dir([subjrootfolder '\*.mat']);
    rootfilenames = cell(1,length(rootfiles));
    rootfilenamescmp = cell(1,length(rootfiles));
    for ifile = 1:length(rootfiles)
        rootfilenames{ifile} = rootfiles(ifile).name;
        
        % Reverse string and remove '.mat', for comparison with strncmpi
        rootfilenamescmp{ifile} = rootfiles(ifile).name(end-4:-1:1); 
    end
    clear rootfiles;
    
    % Check if storage filenames specified by Trials field exist. 
    % If so, find the filename, if not, create new filename
    selfilenames = cell(1,length(subjtrials));
    for itrial = 1:length(subjtrials)
        
        otherStr = subjtrials{itrial};
        matches = strncmpi(otherStr(end:-1:1),rootfilenamescmp,length(subjtrials{itrial}) );
        
        if any(matches)
            selfilenames{itrial} = rootfilenames{matches};
        else
            bckslsh = strfind(subjrootfolder,'\');
            if isempty(bckslsh)
                bckslsh = 0;
            end
            selfilenames{itrial} = [subjrootfolder(bckslsh(end)+1:end) '_' subjtrials{itrial}];
        end
    end
    
    % Obtain inputs, if any
    if ~isempty(get(htabmod,'children'))

%         if any(isubj ~= subjExc)||isempty(subjExc) % From the module fields
        if ~any(isubj == subjExc)||isempty(subjExc) % From the module fields
            
            % This is where things get nasty : dynamic variable creation
            hfields = findobj(htabmod,'-regexp','tag','edt');
            addInputStr = cell(1,length(hfields));
            for ifield = 1:length(hfields) 
                tag = get(hfields(ifield),'tag');
                unders = strfind(tag,'_');
                
                fieldstr = get(hfields(ifield),'string');
                if isempty(fieldstr)
                    fieldstr = '[]';
                end
                
                % Convert to single string (get rid of multiple rows)
                if ~any( size(fieldstr) == 1)
                    fieldstr = regexprep(strjoin(cellstr(fieldstr)),'\.\.\.','');
                end

                eval(['modInput' tag(unders(end)+1:end) ' = ' fieldstr ';']);
                addInputStr{str2double(tag(unders(end)+1:end))} = [',modInput' tag(unders(end)+1:end)];
            end
            addInputStr = strjoin(addInputStr); % Convert from cell to single string

        else % From the exception fields

            % This is where things get nasty : dynamic variable creation
            htabexc = findobj(htabgexc,'tag',['tabExc_' modidx '_' num2str(isubj)]);
            hfields = findobj(htabexc,'-regexp','tag','edt');
            addInputStr = cell(1,length(hfields));
            for ifield = 1:length(hfields)
                tag = get(hfields(ifield),'tag');
                unders = strfind(tag,'_');
                
                fieldstr = get(hfields(ifield),'string');
                if isempty(fieldstr)
                    fieldstr = '[]';
                end
                
                % Convert to single string (get rid of multiple rows)
                if ~any( size(fieldstr) == 1)
                    fieldstr = regexprep(strjoin(cellstr(fieldstr)),'\.\.\.','');
                end
                
                eval(['modInput' tag(unders(end)+1:end) ' = ' fieldstr  ';']);
                addInputStr{str2double(tag(unders(end)+1:end))} = [',modInput' tag(unders(end)+1:end)];
            end
            addInputStr = strjoin(addInputStr); % Convert from cell to single string
            
        end
        
    else
        addInputStr = '';
    end

    % Update subject progress bar
    set( findobj(h1338,'tag','txtProgBarS') , 'position' , [0 0.975 find(isubj==subjects)/length(subjects) 0.025] );
    set( findobj(h1338,'tag','txtSubjProg') , 'string',['Current subject: ' num2str(isubj)])
    pause(0.25); % Give the UI some time to update
    
    % Trial loop
    for itrial = 1:length(subjtrials)
        
        % Load or create trial (subject field is always updated)
        if exist( [subjrootfolder '\' selfilenames{itrial}],'file')
            load([subjrootfolder '\' selfilenames{itrial}]); % load Datastr structure
        end
        Datastr.Info.SubjRoot = subjrootfolder;
        Datastr.Info.Trial = subjtrials{itrial};
        
        % Get other info from 'other' field
        % NOTE: you're not allowed to assign multiple variables in a single line of text in the otherField
        otherStr = cellstr(get(findobj(tabgSubj,'tag',['edtOther' num2str(isubj)]),'string'));
        if any(~cellfun(@isempty,strfind(otherStr,'...')))
            idxCell = 1;
            while any(~cellfun(@isempty,strfind(otherStr,'...')))
                if isempty(strfind(otherStr{idxCell},'...'))
                    idxCell = idxCell + 1;
                else
                    otherStr{idxCell} = regexprep(strjoin(otherStr(idxCell:idxCell+1)),'\.\.\.','');
                    otherStr(idxCell+1) = [];
                end
            end
        end
        
        % Dynamically create and store the variables in other string
        % NOTE: this field is not checked in advance. We'll just skip if an error occurs.
        for idxCell = 1:size(otherStr,1)
            try
                foo = otherStr{idxCell};
                idxIs = strfind(foo,'=');
                eval( [ 'Datastr.Info.(''' regexprep(foo(1:idxIs-1),' ','') ''') = ' regexprep(foo(idxIs+1:end),' ','') ';' ] );
            catch
                warning('gBMP_modrun:other',['Error(s) in evaluationg Other field of subject ' num2str(isubj) '. Skipping.']);
            end
        end
        
        % Update trial progress bar
        set( findobj(h1338,'tag','txtProgBarT') , 'position' , [0 0.95 itrial/length(subjtrials) 0.025] );
        set( findobj(h1338,'tag','txtTrialProg') , 'string',['Current trial: ' subjtrials{itrial}])
        pause(0.25); % Give the UI some time to update
        
        if ~keepLegacy
        
            % Module function call
            % (There is no check on whether the structure contains the fields required for running the module)
            eval( ['Datastr = ' get(htabmod,'title') '(Datastr' addInputStr ');'] );
            
            % Save
            if ~isempty(Datastr)
                save([subjrootfolder '\' selfilenames{itrial}],'Datastr');
            end
        else
            
            % Create copy
            DatastrLeg = Datastr;
            
            % Module function call
            % (There is no check on whether the structure contains the fields required for running the module)
            eval( ['Datastr = ' get(htabmod,'title') '(Datastr' addInputStr ');'] );
            
            % Save
            if ~isempty(Datastr)
                save([subjrootfolder '\' selfilenames{itrial}],'Datastr');

                foo = clock;
                save([subjrootfolder '\' selfilenames{itrial} num2str(foo(end-2)) num2str(foo(end-1))],'DatastrLeg'); % Legacy file
            end
            
        end
        clear Datastr DatastrLeg foo
        
    end
    
end

%% Make all fields visible again

set(findobj(h1338,'tag','tabgSubj'),'visible','on');
set(findobj(h1338,'tag','tabgMod'),'visible','on');
set(findobj(h1338,'-regexp','tag','tabgExc'),'visible','on');
set(findobj(h1338,'tag','psh>'),'enable','on');
gBMP_modenablexcpt(hObject,eventdata,h1338);

set(findobj(h1338,'tag','uifile'),'enable','on');
set(findobj(h1338,'tag','uiproc'),'enable','on');
set(findobj(h1338,'tag','uihelp'),'enable','on');

% Update progress bars
set( findobj(h1338,'tag','txtProgBarS') , 'position' , [0 0.975 0 0.025] );
set( findobj(h1338,'tag','txtSubjProg') , 'string','Current subject: 0');
set( findobj(h1338,'tag','txtProgBarT') , 'position' , [0 0.95 0 0.025] );
set( findobj(h1338,'tag','txtTrialProg') , 'string','Current subject: 0');

end