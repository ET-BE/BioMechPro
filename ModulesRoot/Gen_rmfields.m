function [Datastr] = Gen_rmfields(Datastr,field2rm)
% gBMPDynUI field2rm=1;
%
% INPUT)
% - Datastr: the data structure with any number of fields
% 
% - field2rm: string or cell of strings, specifying the field(s) in Datastr
% that should be removed
% 
% - Datastr: the data structure, without the fields specified by field2rm


%% Check 

if ischar(field2rm)
    field2rm = {field2rm};
end

%% Find the fields to remove

% Find
for fld = 1:numel(field2rm)
    
    labelField = field2rm{fld};
    
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
            warning('Gen_butter:noField',['Cannot find field ' labelField '. Skipping.']);
        end

    else
        myField = layer1(isname1);
    end

    if numel(myField) == 1
        Datastr = rmfield(Datastr,myField{1});
    elseif numel(myField) == 2
        Datastr.(myField{1}) = rmfield(Datastr.(myField{1}),myField{2});
    end
    
end



end