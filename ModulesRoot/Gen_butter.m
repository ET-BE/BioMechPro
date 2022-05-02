function [Datastr] = Gen_butter(Datastr,field2filt,dim2filt,fOrd,fCut,fType,fCall)
% gBMPDynUI field2filt=1; dim2filt=1; fOrd=1; fCut=1; fType=1; fCall=1;
% 
% Filter data in fields along specified dimension using filtfilt or filter
% commands.
% 
% INPUT)
% - Datastr: structure, with any numerical field that can be filtered
% 
% - field2filt: string or cell of strings, specifying the names of the
% fields to filter.
% 
% - dim2filt: integer or vector of integers, specifying the dimensions
% along which the filter operation has to occur. dim2filt must have the 
% same number of elements as field2filt
% 
% - fOrd: integer or vector of integers, specifying the filter order for
% each field specified by filed2filt. fOrd must have the same number of
% elements as field2filt.
% 
% - fCut: integer, vector of integers, or cell of integers and/or vectors, 
% specifying the cutoff frequency (or frequencies for band filters) 
% normalized to the Nyquist frequency, for each field specified by 
% field2filt. 
% If fCut is a vector, it must have the same number of elements as
% field2filt. This means that vector input can only be used for low or high
% pass filters, as band filters require 2 cutoff frequencies. When using
% band filters, use cell notation.
% If fCut is a cell, it must have the same number of cell elements as
% field2filt.
% 
% - fType: string or cell of strings, with one of 'high', 'stop', 'low', or
% 'bandpass', to specify the type of filter.
% 
% - fCall: string or cell of strings, with one of 'filtfilt' or 'filter',
% specyfing the filter operation. If left empty or filled out incorrectly, 
% it will default to 'filtfilt'.
% 
% OUTPUT)
% - Datastr: structure, with filtered fields specified by field2filt.
% 
% NOTES)
% This function doesn't look for fields that are deeper than the 2nd level
% (With the main structure being level 0).
% 
% EXAMPLE)
% Suppose structure Datastr.Marker.MarkerData and Datastr.Force.ForceData,
% both with numerical data, then:
% field2filt = {'MarkerData','ForceData'};
% dim2filt = [1 1];
% fOrd = [1 2];
% fCut = {20/100,[40/100 60/100]}
% fType = {'low','stop'};
% fCall = {'filtfilt','filter'}
% Gen_butter(Datastr,field2filt,dim2filt,fOrd,fCut,fType,fCall) gives:
% 1) Datastr.Marker.MarkerData, with data filtered along the first dimension 
% with a zero phase 2nd order 20 Hz low pass (note: zero phase doubles
% order: see doc filtfilt). 
% 2) Datastr.Force.ForceData, with data filtered along the first dimension
% with a 2nd order 40-60 Hz bandstop digital filter along the first dimension.
% 

% TODO)
% When you encounter a field that's incorrect or cannot be filtered (e.g.
% warning statement + return), then all remaining fields are skipped,
% rather than just skipping the incorrect one.

%% Do a check on the inputs

if ischar(field2filt)
    field2filt = {field2filt};
end

% Check dim2filt
if numel(dim2filt)~=numel(field2filt)
    warning('Number of elements in dim2filt ~= number of elements in field2filt. Skipping');
    return;
end

% Check fOrd
if numel(fOrd)~=numel(field2filt)
    warning('Number of elements in fOrd ~= number of elements in field2filt. Skipping');
    return;
end

% Check fCut
if isnumeric(fCut)
    fCut = num2cell(fCut);
elseif numel(fCut)~=numel(field2filt)
    warning('Number of elements in fCut ~= number of elements in field2filt. Skipping');
    return;
end

% Check fType
if ischar(fType)
    fType = {fType};
end
if numel(fType)~=numel(field2filt)
    warning('Number of elements in fType ~= number of elements in field2filt. Skipping');
    return;
end

% Check fCall
if isempty(fCall)
    fCall = {'filtfilt'}; % Default
elseif ischar(fCall)
    fCall = {fCall};
end
if numel(fCall)~=numel(field2filt)
    warning('Number of elements in fCall ~= number of elements in field2filt. Skipping');
    return;
end

%% Find the label fields to operate on

for fld = 1:length(field2filt)
    
    labelField = field2filt{fld};
    
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

    
    
    %% Filter the field
    
    if ~ isempty(myField)
        
        % Do some checks
        if numel(myField) == 1
            dat2filt = Datastr.(myField{1});
        elseif numel(myField) == 2
            dat2filt = Datastr.(myField{1}).(myField{2});
        end
        if ~isnumeric(dat2filt)
            warning('Gen_diff:nonumeric',['Data in field ' field2filt{fld} ' is not numeric. Skipping.']);
            return; % TODO : this also skips all other fields in field2diff. Should go back to the start of the loop
        elseif ndims(dat2filt) < dim2filt(fld)
            warning('Gen_diff:dims',['Number of dimensions of data in ' field2filt{fld} ' < specified dimension. Skipping.']);
            return;
        end

        % Create filter
        [b,a] = butter(fOrd(fld),fCut{fld},fType{fld});

        % Filter data
        if strcmpi(fCall{fld},'filter')
            dat2filt = filter(b,a,dat2filt,[],dim2filt(fld));
        else

            % Permute the data if required (only for filtfilt)
            if dim2filt(fld)~=1
                dimo = setdiff(dim2filt(fld),ndims(dat2filt));
                dat2filt = permute(dat2filt,[dim2filt(fld) dimo]);
            end

            dat2filt = filtfilt(b,a,dat2filt);

            % Return data to its original dimensions, if required
            if dim2filt(fld)~=1
                [~,perm] = sort([dim2filt(fld) dimo],'ascend');
                dat2filt = permute(dat2filt,perm);
            end

        end

        % Store data
        if numel(myField) == 1
            Datastr.(myField{1}) = dat2filt;
        elseif numel(myField) == 2
            Datastr.(myField{1}).(myField{2}) = dat2filt;
        end
    
    end
    
end


end