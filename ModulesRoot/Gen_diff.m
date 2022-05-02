function [Datastr] = Gen_diff(Datastr,field2diff,ndiff,dim2diff,fs)
% gBMPDynUI field2diff=1; ndiff=1; dim2diff=1; fs=1;
% 
% Differentiate data in specified fields along specified dimensions
% 
% INPUT)
% - Datastr: structure, with any numerical field that can be differentiated
% 
% - field2diff: string or cell of strings, specifying the names of the
% fields to be differentiated.
% 
% - dim2diff: integer or vector of integers, specifying the dimensions
% along which the differentiation is to occur. Dim2diff must have the same
% size as field2diff.
% 
% - fs: integer or vector of integers, specifying the sample frequency (Hz)
% with which the differentiated data is multiplied
% 
% - ndifs: integer or vector of integers, specifying the number of times
% the data has to be differentiated. 
% 
% OUTPUT)
% - Datastr: structure, with added fields having names specified by
% field2diff, and an additional 'D' for each differential step.
% 
% NOTES)
% The differentiated data will be padded by zeros at the start, to ensure
% the same lenght as the original data.
% 
% This function doesn't look for fields that are deeper than the 2nd level
% (With the main structure being level 0).
% 
% EXAMPLE)
% Suppose structure Datastr.Marker.MarkerData with numerical data, then:
% Gen_diff(Datastr,'MarkerData',2,1,100)
% Will differentiate the field MarkerData along the first dimension, twice,
% generating the fields:
% Datastr.Marker.MarkerDataD
% Datastr.Marker.MarkerDataDD
% which contain the first and second derivatives respectively.

%% Do a check on the inputs

if ischar(field2diff)
    field2diff = {field2diff};
end

if isempty(ndiff)
    ndiff = ones(size(field2diff));
elseif numel(ndiff)~=numel(field2diff)
    warning('Gen_diff:numel1','Number of elements in ndiff ~= number of elements in field2diff. Skipping.');
    return;
end
if isempty(dim2diff)
    dim2diff = ones(size(field2diff));
elseif numel(dim2diff)~=numel(field2diff)
    warning('Gen_diff:numel2','Number of elements in dim2diff ~= number of elements in field2diff. Skipping.');
    return;
end
if isempty(fs)
    fs = ones(size(field2diff));
elseif numel(fs)~=numel(field2diff)
    warning('Gen_diff:numel','Number of elements in fs ~= number of elements in field2diff');
    return;
end

% For fieldname differentiated data
addD = 'D';

%% Find the label fields to operate on
for fld = 1:length(field2diff)
    
    labelField = field2diff{fld};
    
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
            warning('Gen_diff:noField',['Cannot find field ' labelField '. Skipping.']);
        end

    else
        myField = layer1(isname1);
    end
    
    
    %% Differentiate the field (if numeric), and dynamically add a new one
    
    if ~isempty(myField)
    
        % Do some checks
        if numel(myField) == 1
            dat2diff = Datastr.(myField{1});
        elseif numel(myField) == 2
            dat2diff = Datastr.(myField{1}).(myField{2});
        end
        if ~isnumeric(dat2diff)
            warning('Gen_diff:nonumeric',['Data in field ' field2diff{fld} ' is not numeric. Skipping.']);
            return;
        elseif ndims(dat2diff) < dim2diff(fld)
            warning('Gen_diff:dims',['Number of dimensions of data in ' field2diff{fld} ' < specified dimension. Skipping.']);
            return;
        end

        diffdat = zeros(size(dat2diff));
        for nd = 1:ndiff(fld)

            foo = diff(dat2diff,nd,dim2diff(fld))*(fs(fld).^nd);
            dimo = setdiff(1:ndims(foo),dim2diff(fld)); % All other dimensions

            % Selection tool for all other dimensions
            if ~ isempty(dimo)
                sel = cell(1,length(dimo));
                for idx = 1:length(dimo)
                    sel{idx} = ':';
                end
            else
                sel = {':'};
            end

            % Swap dimensions to padd user specified dimension, then swap back
            foo = permute(foo,[dim2diff(fld) dimo]);
            diffdat = permute(diffdat,[dim2diff(fld) dimo]);
            diffdat(1+nd:end,sel{:}) = foo;
            [~,perm] = sort([dim2diff(fld) dimo],'ascend');
            diffdat = permute(diffdat,perm);

            % Store data
            if numel(myField) == 1
                Datastr.([myField{1} repmat(addD,[1 nd])]) = diffdat;
                Datastr = orderfields(Datastr);
            elseif numel(myField) == 2
                Datastr.(myField{1}).([myField{2} repmat(addD,[1 nd])]) = diffdat;
                Datastr.(myField{1}) = orderfields(Datastr.(myField{1}));
            end

        end

    end
    
end

end