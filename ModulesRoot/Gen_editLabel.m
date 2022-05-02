function [Datastr] = Gen_editLabel(Datastr,labelField,cellStrOld,cellStrNew,keepInput)
% gBMPDynUI labelField=1; cellStrOld=1; cellStrNew=1; keepInput=1
% 
% Each element of the cell with name labelField present in structure Datastr
% must contain a string. Each instance of cellStrOld within each of those
% strings will be replaced by each instance of cellStrNew.
% For example, if within input Datastr you have:
% Datastr.Marker.MarkerDataLabel = {'aA','aB','bA','bB'}
% Then for inputs:
% labelField = 'MarkerDataLabel';
% cellStrOld = {'a','b'}
% cellStrNew = {'b','a'}
% You will get:
% Datastr.Marker.MarkerDataLabel = {'bA','bB','aA','aB'}
% 
% INPUT)
% 
% - Datastr: structure, containing a cell you want to edit
% 
% - labelField: string, fieldname containing the cell you want to edit.
% Example:
% 'MarkerDataLabel'
% To edit the field present in Datastr.Marker.MarkerDataLabel
% 
% - cellStrOld: cell, containing (part of) strings in field with name 
% labelField that need to be replaced.
% 
% -cellStrNew: cell, containing strings that will replace those instances
% of cellStrOld.
% 
% -keepInput: boolean, if TRUE, the input will be stored in the structure
% such that it can be re-used, for example in probe reconstruction (where
% you also have to swap the labels in the probe trials)
% 
% OUTPUT)
% Datastr: structure, with adjusted labelcell
% 
% NOTES)
% cellStrOld and cellStrNew must have the same number of elements
% 
% cellStrOld may not contain two entries that are the same
% 
% All inputs are case sensitive
% 
% This function doesn't look for fields that are deeper than the 2nd level
% (With the main structure being level 0).

%% Check input

if isempty(cellStrOld) || isempty(cellStrNew) || isempty(labelField)
    return; % Skip without warning
end

if ischar(cellStrNew)
    cellStrNew = {cellStrNew};
end
if ischar(cellStrOld)
    cellStrOld = {cellStrOld};
end

if length(cellStrOld) ~= length(cellStrNew)
    warning('Gen_swapLabel:size','cellStrOld and cellStrNew must have same number of elements.');
    return;
end

%% Find the label field to operate on

layer1 = fieldnames(Datastr);

isname1 = strcmp(layer1, labelField); % Mind that strcmp is used here (case sensitive)
if ~any(isname1)
    
    myField = {};
    for ifld = 1:length(layer1)

        layer2 = fieldnames( Datastr.(layer1{ifld}) );
        isname2 = strcmp(layer2, labelField);
        if any(isname2)
            myField = {layer1{ifld} , layer2{isname2}};
            break
        end
        
    end
    if isempty(myField)
        warning('Gen_swapLabel:noField',['Cannot find field ' labelField '. Skipping.']);
        return;
    end
    
else
    myField = {layer1(isname1)};
end


%% Change the labelcell

% Get labelcell
if numel(myField) == 1
    labelcell = Datastr.(myField{1});
elseif numel(myField) == 2
    labelcell = Datastr.(myField{1}).(myField{2});
end

if ~iscell(labelcell)
    warning('Gen_swapLabel:noCell','field labelField is not a cell. Skipping.')
    return;
end


% Edit the cell
newParts = cell(numel(cellStrOld)+1,numel(labelcell)); % +1 for the original cell
newParts(1,:) = labelcell;
for istr = 1:numel(cellStrOld)
    
    newParts(istr+1,:) = regexprep(labelcell,cellStrOld{istr},cellStrNew{istr});
    
end

% The least occuring string in the columns is the one you want
labelcellnew = cell(size(labelcell));
for ilbl = 1:size(newParts,2)
    
    col = newParts(:,ilbl);
    [uni,~,idx] = unique(col,'stable');
    [~,minidx] = min( histcounts(idx,length(uni)) );
    
    labelcellnew(ilbl) = uni(minidx);
    
end

% Assign the new labelcell
if numel(myField) == 1
    Datastr.(myField{1}) = labelcellnew;
elseif numel(myField) == 2
    Datastr.(myField{1}).(myField{2}) = labelcellnew;
end


% Store the input for re-use in other modules
if keepInput
    
    Datastr.Info.editLabel2 = labelField;
    Datastr.Info.editLabel3 = cellStrOld;
    Datastr.Info.editLabel4 = cellStrNew;
    
end

end