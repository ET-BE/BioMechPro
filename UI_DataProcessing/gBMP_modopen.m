function [] = gBMP_modopen(hObject,eventdata,h1338)
% Open module and dynamically create user input fields in the module tab
% when field info is specified in the module function.
% 
% Fields are spaced over normalized tab height range 0.55-1
% The default height of a single field normalized to tab height is 0.068
% The default space between fields and tab boundaries is 0.0097
% So the size of a single field is 0.068*nlines + 0.0097
% If the fields don't fit in a column, shrink the field size

space2distr = 1;
fheight = 0.068;
fspace = 0.0097;

% Find handle to module tabgroup
htabgmod = findobj(h1338,'tag','tabgMod');

% Find the program's root folder on the Matlab search path
allpaths = path;
[tok,remtok] = strtok(allpaths,';');
if isempty(strfind(tok,'BioMechPro\ModulesRoot'))
    while ~isempty(remtok)
        [tok,remtok] = strtok(remtok,';');
        if ~isempty(strfind(tok,'BioMechPro\ModulesRoot'))
            break
        end
    end
end

% Let user select module(s) (default in module folder if possible)
if ~isempty(tok)
    [filenames,pathname] = uigetfile('.m','Select processing module(s)',...
        [tok '\'],...
        'multiselect','on');
else
    [filenames,pathname] = uigetfile('.m','Select processing module(s)',...
        'multiselect','on');
end

% Return if selection is cancelled
if isnumeric(filenames)
    return;
end

% If single selection, convert to cell to keep format same as in multi-select
% (pathnames is always a string: cannot multi-select in multiple folders)
if ~iscell(filenames) % Single select
    foo{1} = filenames;
    filenames = foo;
end

% Check if any module is already opened
if strcmpi(get(htabgmod,'visible'),'on')

    % Check if any of the selected modules is already open. If so, ignore.
    for ifile = 1:length(filenames)
        filename = filenames{ifile};
         
        if ~isempty( findobj(h1338,'parent',htabgmod,'-and','title',filename(1:end-2)) ) % Ignore if already opened, -2 to get rid of the .m
            filenames{ifile} = []; % Mark for removal
        end
    end
    
end

% Remove empty cells and return if nothing is left
filenames(cellfun('isempty',filenames)) = []; 
if isempty(filenames)
    return;
end

%% Dynamically create UI tabs and UI fields in those tabs
for ifile = 1:length(filenames)
    
    filename = filenames{ifile};
    
    % Read info from function file for UI edit field creation, if any
    try
        funtxt = fileread([pathname filename]);
        idxstart = strfind(funtxt,'gBMPDynUI') + 9;
        idxend = strfind(funtxt,char(10)); % char(10) is the newline character
        idxend = idxend( find((idxend - idxstart)>0,1,'first') );
        uiinfo = regexprep( funtxt(idxstart:idxend) ,' ','');
        uiinfo = textscan(uiinfo,'%s','Delimiter',';'); % uiinfo is now a cell containing another cell
        uiinfo = uiinfo{1}; 
    catch
%         warning('gBMP_modopen:fopen',['Cannot read ' filename ' for info to create UI fields, if any.']);
        uiinfo = '';
    end

    % Determine field distribution inside tab
    nlines = zeros(1,length(uiinfo));
    for ifield = 1:length(uiinfo) 
        finfo = uiinfo{ifield};
        
        idxsplit = strfind(finfo,'=');
        fieldname{ifield} = finfo(1:idxsplit-1);
        nlines(ifield) = str2double(finfo(idxsplit+1:end));
    end

    if ~isempty(nlines)
        % Determine field order inside tab: always place largest first
        [nlinessort,permord] = sort(nlines,'descend');
        fieldcol = zeros(size(nlinessort));
        colsiz1 = 0; colsiz2 = 0;
        for idxsel = 1:length(nlinessort)

            if colsiz1 <= colsiz2   % Place field in first column
                colsiz1 = colsiz1 + nlinessort(idxsel)*fheight + fspace;
                fieldcol(idxsel) = 1;
            else                    % place field in second column
                colsiz2 = colsiz2 + nlinessort(idxsel)*fheight + fspace;
                fieldcol(idxsel) = 2;
            end

        end

        % Scale field size if needed
        if colsiz1 > space2distr
            nlinessort(fieldcol==1) = nlinessort(fieldcol==1) * (space2distr/colsiz1);
        end
        if colsiz2 > space2distr
            nlinessort(fieldcol==2) = nlinessort(fieldcol==2) * (space2distr/colsiz2);
        end
        
    end
    
    % Get existing tabs
    modtabs = get(htabgmod,'children');
    
    % Find an index to assign to the module tab
    if isempty(modtabs)
        modidx = 1;
    else
        % Find the lowest positive integer not among existing tab indices
        modidx = zeros(1,length(modtabs));
        for ichild = 1:length(modtabs)
            foo = get(modtabs(ichild),'tag');
            unders = strfind(foo,'_');
            modidx(ichild) = str2double(foo(unders(end)+1:end));
        end
        modidx = min(setdiff(1:length(modtabs)+1,modidx));
    end

    % Create module tab
    tabMod = uitab(htabgmod,...
        'backgroundcolor',[0.9 0.9 0.95],...
        'foregroundcolor',[0 0 0],...
        'title',filename(1:end-2),...
        'tag',['tabMod_' num2str(modidx)]);
    
    % Create fields inside module tab
    for ifield = 1:length(uiinfo)
        
        if sum(fieldcol(1:ifield)==fieldcol(ifield)) == 1 % If first in column
            vertipostxt = 1 - fheight - fspace*sum(fieldcol(1:ifield)==fieldcol(ifield)) ;
        else
            vertipostxt = 1 - fheight - fspace*sum(fieldcol(1:ifield)==fieldcol(ifield)) - fheight*sum(nlinessort(fieldcol(1:ifield-1)==fieldcol(ifield)));
        end       
        vertiposedt = 1 - fspace*sum(fieldcol(1:ifield)==fieldcol(ifield)) - fheight*sum(nlinessort(fieldcol(1:ifield)==fieldcol(ifield)));
        horizpos = 0.01+0.5*(fieldcol(ifield)-1);
        
        % Text object (fieldname)
        uicontrol(tabMod,...
        'backgroundcolor',[0.9 0.9 0.95],...
        'foregroundcolor',[0 0 0],...
        'fontsize',10,...
        'fontweight','bold',...
        'horizontalalignment','left',...
        'string',fieldname(permord(ifield)),...
        'style','text',...
        'tag',['txtMod_' num2str(modidx) '_' num2str(permord(ifield))],...
        'units','normalized',... % Be sure to set units before position
        'position',[horizpos vertipostxt 0.1 fheight]);
    
        % Edit object
        uicontrol(tabMod,...
        'backgroundcolor',[1 1 1],...
        'foregroundcolor',[0 0 0],...
        'fontsize',10,...
        'horizontalalignment','left',...
        'max',nlines(permord(ifield)),...
        'string','',...
        'style','edit',...
        'tag',['edtMod_' num2str(modidx) '_' num2str(permord(ifield))],...
        'units','normalized',... % Be sure to set units before position
        'position',[horizpos+0.1 vertiposedt 0.34 fheight*nlinessort(ifield)]);
        
    end
    
    % Exception tabgroup (on main figure!)
    uitabgroup(h1338,...
        'position',[0 0 1 0.35],...  % Normalized on primary module tab
        'tag',['tabgExc_' num2str(modidx)],...
        'visible','off');
end

% Make module tab visible and enable run button
if strcmpi(get(htabgmod,'visible'),'off')
    set(htabgmod,'visible','on');
    set(findobj(h1338,'tag','psh>'),'enable','on');
end

% Switch to last created module tab and run the 'tabchanged' function
set(htabgmod,'selectedtab',tabMod);
gBMP_modenablexcpt(hObject,eventdata,h1338);

end