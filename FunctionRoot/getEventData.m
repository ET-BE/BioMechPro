function seldata = getEventData(data,event,dims,dimr,varargin)
%% getEventData
% Get (time series) data at events specified by an event vector or matrix
% 
% INPUT)
% -data : vector or matrix, containing (time series) data. 
% -event : vector or matrix containing only zeros (or false), except a value 1 (or true)
% at instances of an event at which data in "data" must be selected. 
% The dimensions of "event" must correspond with the length of "data" in 
% dimensions "dims" and "dimr"
% -dims : scalar, the dimension in data in which the samples (over time) are located
% -dimr : scalar, the dimension in data in which repetitions (of time series) are located
% If data is a vector without repetitions, specify dimr = 1
% (The intended use of this function is for multiple time series in a single matrix).
% 
% -nevt : scalar, number of events to select from "event". Must be less than
% or equal to the number of events in each repetition of "event".
% Specify when "event" does not contain an equal amount of events for each repetition.
% Selection occurs from first to last.
% 
% OUTPUT)
% seldata : vector or matrix, with data at the selected events.
% The size of "seldata" is equal to the size of input "data", 
% except for dimension "dims" which is of length "nevent".
% 
% NOTES)
% 
% 
% EXAMPLE)
% Some test data
% data = rand(1,2,3,4,5);
% event = zeros(3,5);
% event(1,3) = 1; event(2,4:5) = 1; event(3,2) = 1;
% dims = 5;
% dimr = 3;
% nevt = 1;
% 
% foo = getEventData(data,event,dims,dimr,nevt);
% 
% gives you:
% foo = data(:,:,2,:,[4 5])
% 
% The first and last row of event are not used, as they contain fewer events
% than specified by nevt.
% 
% If instead you would use nevt = 1, you would get a single matrix foo containing:
% data(:,:,1,:,3)
% data(:,:,2,:,4)
% data(:,:,3,:,2)
% 
% Which you can check from:
% foo(:,:,1,:) - data(:,:,1,:,3)
% foo(:,:,2,:) - data(:,:,2,:,4)
% foo(:,:,3,:) - data(:,:,3,:,2)
% 
% which will give you only zeros

%% Check input
ndim = ndims(data);
dsiz = size(data);

% Make event logical
event = event ~= 0;

% Varargin
if nargin == 4
    nevt = [];
elseif nargin > 4
    nevt = varargin{1};
elseif nargin > 5
    error('getEventData: too many input arguments');
end

% General dimensions
if (dims > ndim) || (dimr > ndim)
    error('getEventData: dimension dims and/or dimr are/is too high');
end

% Event dimensions vs Data dimensions
if (size(event,2) == dsiz(dims))&&(size(event,1) == dsiz(dimr))
    event = event';
elseif (size(event,1) ~= dsiz(dims))||(size(event,2) ~= dsiz(dimr))
    error('getEventData: dimensions of event should correspond with those in data, as specified by dims and dimr');
end

% Number of events
if (dimr > 1) && isempty(nevt)
    if ( any(diff(sum(event,1))~=0) ) 
        error('getEventData: parameter event does not contain an equal number of events for all repetitions. Specify parameter nevt.');
    else
        nevt = sum(event(:,1)); % number of events
    end
elseif (dimr > 1) && ( any(nevt > sum(event,1)) )
    % Option 1: generate warning and remove insufficient ones (see below)
%     warning('getEventData: number of events in event input is sometimes less than events than specified by nevt. These will be removed.');
    
    % Option 2: generate error
%     error('getEventData: number of specified events cannot be larger than number of events in event input');
elseif dimr == 1
    if isempty(nevt)
        nevt = sum(event);
    elseif nevt > sum(event)
        error('getEventData: number of specified events cannot be larger than number of events in event input');
    end
end

%% Do stuff

% Find dimensions other than dims or dimr (if any)
dimo = setdiff(1:ndim,[dims dimr]);

% Create a selection tool for all data in dimensions other than dims or dimr
if ~isempty(dimo)
    sel = cell(1,length(dimo));
    for i = 1:length(sel)
        sel{i} = ':';
    end
else
    sel = {':'};
end

% Permute data with dims and dimr in first and second dimension
% (to generalize dimensions, makes writing code easier)
pdata = permute(data,[dims dimr dimo]);

% Remove repetitions with insufficient events (in case of warning)
pdata(:, nevt > sum(event,1) , sel{:}) = [];
event(:, nevt > sum(event,1)) = [];

% Select the first nevt events
eventc = cumsum(event,1);
eventc(~event) = inf;
event = eventc <= nevt;

% WARNING: THIS ONLY WORKS IF EVERY EVENT COLUMN HAS THE SAME NUMBER OF EVENTS!
% Replicate event over the other dimensions, make logical
eventr = repmat(event,[1 1 dsiz(dimo)]);

% Select events in repetition
seldata = reshape( pdata(eventr) , [nevt size(event,2) dsiz(dimo)] );

% Return data to original dimensions
[~,perm] = sort([dims dimr dimo],'ascend');
seldata = permute(seldata,perm);

end