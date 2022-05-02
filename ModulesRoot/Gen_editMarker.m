function [Datastr] = Gen_editMarker(Datastr,dataField,dim,markIdxOld,markIdxNew,keepInput)
% gBMPDynUI dataField=1; dim=1; markIdxOld=1; markIdxNew=1; keepInput=1;
% 
% INPUT)
% 
% - Datastr: structure, containing a cell you want to edit
% 
% - dataField: string, fieldname containing the matrix you want to edit.
% Example:
% 'MarkerData'
% To edit the field present in Datastr.Marker.MarkerData
% 
% - dim: scalar, the dimension in which the swap has to take place
% 
% - markIdxOld: vector, containing indices of the markers in dataField that 
% need to be replaced.
% 
% - markIdxNew: vector, containing indices of the markers in dataField that
% will replace those markers specified by markIdxOld.
% 
% -keepInput: boolean, if TRUE, the input will be stored in the structure
% such that it can be re-used, for example in probe reconstruction (where
% you also have to swap the labels in the probe trials)
% 
% OUTPUT)
% Datastr: structure, with adjusted labelcell
% 
% NOTES)
% 
% The swap always works on the 2nd dimension
%
%

%% Do some checks

if isempty(dataField)||isempty(dim)||isempty(markIdxOld)||isempty(markIdxNew)
    return; % Return without warning
end

if numel(markIdxOld)~=numel(markIdxNew)
    warning('Gen_editMarker:numel','markIdxOld and markIdxNew must have same number of elements. Skipping.')
    return;
end

%% Find the data field to operate on

layer1 = fieldnames(Datastr);

isname1 = strcmp(layer1, dataField); % Mind that strcmp is used here (case sensitive)
if ~any(isname1)
    
    myField = {};
    for ifld = 1:length(layer1)

        layer2 = fieldnames( Datastr.(layer1{ifld}) );
        isname2 = strcmp(layer2, dataField);
        if any(isname2)
            myField = {layer1{ifld} , layer2{isname2}};
            break
        end
        
    end
    if isempty(myField)
        warning('Gen_swapLabel:noField',['Cannot find field ' dataField '. Skipping.']);
        return;
    end
    
else
    myField = {layer1(isname1)};
end


%% Swap the markers

% Get data
if numel(myField) == 1
    myData = Datastr.(myField{1});
elseif numel(myField) == 2
    myData = Datastr.(myField{1}).(myField{2});
end

% Some more checks
if ~isnumeric(myData)
    warning('Gen_editMarker:noData','Data in field dataField is not numerical. Skipping.')
    return;
end
if (markIdxOld > size(myData,dim)) || (markIdxNew > size(myData,dim))
    warning('Get_editMarker:dimsiz','markIdxOld or markIdxNew exceeds matrix dimension. Skipping.');
    return;
end

% Create selection variable
sel = repmat({':'},[1 length(size(myData))-1]);

% Other dimensions
dimo = setdiff(1:length(size(myData)),dim);

% Permute data to operate on first dim
myData = permute(myData,[dim dimo]);

% Swap data
myData(markIdxOld,sel{:}) = myData(markIdxNew,sel{:});

% Permute data back
[~,perm] = sort([dim dimo],'ascend');
myData = permute(myData,perm);

% Store data in structure
if numel(myField) == 1
    Datastr.(myField{1}) = myData;
elseif numel(myField) == 2
    Datastr.(myField{1}).(myField{2}) = myData;
end

% Store input
if keepInput
    
    Datastr.Info.editMarker2 = dataField;
    Datastr.Info.editMarker3 = dim;
    Datastr.Info.editMarker4 = markIdxOld;
    Datastr.Info.editMarker5 = markIdxNew;
    
end


end