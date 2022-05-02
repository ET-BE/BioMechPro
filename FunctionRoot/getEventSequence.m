function [redata] = getEventSequence(data,event,dims,dimr,nsmpl,varargin)
%% getEventSequence
% Select and resample input data between two events using interpolation
% 
% INPUT)
% -data : array, containing (time series) data. 
% -event : array consisting of mainly zeros, only containing 
% non-zero values at indices corresponding to a certain event.
% Size must be equal to the lenght of "data" in dimensions "dims" and "dimr".
% If "evts" is not specified, data is interpolated between the first and last
% event for each repetition.
% -dims : scalar, the dimension of data containing the (time series) data.
% -dimr : scalar, the dimension of data in which repetitions (of time series) are located
% If data is a vector without repetitions, specify dimr = 1
% (The intended use of this function is for multiple time series in a single matrix).
% -nsmpl : number of samples to resample
% -evts : vector containing the start and stop number of the events between which to interpolate
% If not specified, interpolation occurs between the first and last event.
% 
% OUTPUT)
% Resampled data between the given events
% 
% NOTES)
% 
% EXAMPLE)
% Suppose we have 12 repetitions of some time series data, each consisting 
% of 500 time samples, in 3 dimensions.
% 
% size(data);
% > [500,3,12]
% 
% Suppose that we have 5 events occuring in the time series data.
% In each repetition these events occur at different time instances.
% This is specified in an event matrix. It must be that:
% 
% size(event);
% >[500,12] (or 12,500)
% 
% We now want to resample the data between events 1 and 2 to 100 samples.
% (i.e. between the indices of the first and second non-zero value in each repetition of event)
% Then:
% 
% redata = getEventSequence(data,event,1,3,100,[1 2])


%% Check input
ndim = ndims(data);
dsiz = size(data);

% Make event logical
event = event ~= 0;

if nargin < 6
    evts = [];
elseif nargin == 6
    evts = sort(varargin{1} , 'ascend');
else
    error('getEventSequence: too many input arguments');
end

% Event dimensions vs Data dimensions
if (size(event,2) == dsiz(dims))&&(size(event,1) == dsiz(dimr))
    event = event';
elseif (size(event,1) ~= dsiz(dims))||(size(event,2) ~= dsiz(dimr))
    error('getEventData: dimensions of event should correspond with those in data, as specified by dims and dimr');
end

% Number of events
if any(evts(1) > sum(event,1)) || any(evts(end) > sum(event,1))
%     warning('getEventSequence: selected events are not (always) available in event variable. These will be removed.');
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
if isempty(evts)
    pdata(:, sum(event,1) < 2 , sel{:}) = [];
    event(:, sum(event,1) < 2) = [];
    % Maybe add warning
else
    pdata(:, sum(event,1) < evts(1) , sel{:}) = [];
    event(:, sum(event,1) < evts(1)) = [];
    
    pdata(:, sum(event,1) < evts(end) , sel{:}) = [];
    event(:, sum(event,1) < evts(end)) = [];
end

% Remove event other than the first and last selected ones
if ~isempty(evts)
    eventc = cumsum(event);
    eventc(~event) = inf;
    event = logical((eventc == evts(1)) + (eventc == evts(end)));
end

% Create storage container
redata = zeros([nsmpl,size(event,2),dsiz(dimo)]);

% Resample every sequence
for irep = 1:size(event,2)
    idxstart = find(event(:,irep),1,'first');
    idxend = find(event(:,irep),1,'last');
    xi = linspace(idxstart,idxend,nsmpl);
   
    redata(1:nsmpl,irep,sel{:}) = interp1(idxstart:idxend,pdata(idxstart:idxend,irep,sel{:}),xi);
end

% Return data to original dimensions
[~,perm] = sort([dims dimr dimo],'ascend');
redata = permute(redata,perm);

end