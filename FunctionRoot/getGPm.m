function [varargout] = getGPm(varargin)
% INPUT)
% Input can be either:
% 
% feetdata : Nx4x3 matrix with data for left and right heel and toe
% The second dimension is assumed as: 
% 1) Left heel
% 2) Left toe
% 3) Right heel
% 4) Right toe
%
% The third dimension is assumed as:
% 1) x direction (right / left)
% 2) y direction (forward / backward : walking direction)
% 3) z direction (up / down)
% 
% or :
% C3DProc : data structure containing at least the fields:...
% feetcell : cell array containing the names of the toe and heel LEDs to
% use. Order is assumed as above (1,2,3,4).
% 
% OUTPUT)
% phasevector : a matrix containing vectors for each phase.
% Non-zero indices indicate the END of each respective phase.
% The value at these indices indicates the amount of samples the phase lasted.
% These vectors are useful when calculating mean / average phase duration
% 
% phasevectorlabel : a labelcell indicating which vector in phasevector corresponds to which gait phase
% NC : no contact
% SSL : single support left
% SSR : single support right
% DS : double support
% DSLIFOR : sub phase, double support left in front of right leg
% DSRIFOL : sub phase, double support right in front of left leg
% 
% phaseidx : structure with the indices of the start and end of each gait phase
% 
% NOTE : requires peakdet.m

% ##############################
% Tuning variable for peakdet.m
% Steps less than this size will not be detected
% TO DO: make this an input parameter in the UI

% peakdetDelta = 0.05; % Meter
peakdetDelta = 0.1; % Meter

% ##############################


%% Check input data

if nargin == 1

    if nargout > 3
        error('getGPm:nout3','Too many output arguments');
    end
    
    feetdata = varargin{1};
    
    if ~isnumeric(feetdata)
        error('getGPm:datatype','Incorrect input type');
    end
    
    if (size(feetdata,2) ~= 4) || (size(feetdata,3) ~= 3)
        error('getGPm:datasiz','incorrect input data size');
    end
    
elseif nargin == 2
    
    if nargout > 1
        error('getGPm:nout1','Too many output arguments');
    end
    
    C3Ddata = varargin{1};
    feetCell = varargin{2};
    
    if ~isstruct(C3Ddata)
        error('getGPm:input1','First of two inputs must be a C3D data structure');
    end
    
    if ~iscell(feetCell) || length(feetCell) ~= 4
        error('getGPm:input2','Second of two inputs must be a 4 element cell with label names of 1) left heel, 2) left toe, 3) right heel, 4) right toe');
    end
    
    % Initiate variables
    labelall = {};
    markerall = [];

    % Get all labels and data
    if isfield(C3Ddata.Marker,'MarkerDataLabel')
        labelall(end+1:end+length(C3Ddata.Marker.MarkerDataLabel)) = C3Ddata.Marker.MarkerDataLabel;
        markerall(:,end+1:end+length(C3Ddata.Marker.MarkerDataLabel),:) = C3Ddata.Marker.MarkerData;
    end
    if isfield(C3Ddata.Marker,'ProbedDataLabel');
        labelall(end+1:end+length(C3Ddata.Marker.ProbedDataLabel)) = C3Ddata.Marker.ProbedDataLabel;
        markerall(:,end+1:end+length(C3Ddata.Marker.ProbedDataLabel),:) = C3Ddata.Marker.ProbedData;
    end

    if isempty(markerall)
        error('getGPm:markerall','No MarkerDataLabel or ProbedDataLabel found in input structure.');
    end
    
    % Get required data
    feetdata = zeros( size(C3Ddata.Marker.ProbedData,1) , 4 , 3 );
    for i = 1:length(feetCell)
        idxmark = strcmpi(feetCell{i},labelall);
        
        if sum(idxmark) == 0
            error('getGPm:feetdata','No data found with name(s) corresponding to those in feetCell');
        end
        
        feetdata(:,i,:) = markerall(:,idxmark,:);
    end
    
else
    error('getGPm:inputs','Incorrect number of inputs');
end


%% Split input data
% 1) Left heel
% 2) Left toe
% 3) Right heel
% 4) Right toe
lhee = squeeze(feetdata(:,1,:));
ltoe = squeeze(feetdata(:,2,:));
rhee = squeeze(feetdata(:,3,:));
rtoe = squeeze(feetdata(:,4,:));

    
%% Find phase onsets and stops

    % Discover events
    [HSL,~] = peakdet(lhee(:,2),peakdetDelta);
    [~,TOL] = peakdet(ltoe(:,2),peakdetDelta);
    [HSR,~] = peakdet(rhee(:,2),peakdetDelta);
    [~,TOR] = peakdet(rtoe(:,2),peakdetDelta);
    % (Note: this doesn't work for backward walking, or "strange" movements)
    
    % Generate error if no events detected
    if any([isempty(HSL) isempty(TOL) isempty(HSR) isempty(TOR)])
%         error(['getGPm:No heel strikes and/or toe-offs detected. Code ' num2str([isempty(HSL) isempty(TOL) isempty(HSR) isempty(TOR)])]);
        warning(['getGPm:No heel strikes and/or toe-offs detected. Code ' num2str([isempty(HSL) isempty(TOL) isempty(HSR) isempty(TOR)]) '. Skipping.']);
        
        if nargin == 1
            varargout{1} = [];
            varargout{2} = [];
            varargout{3} = [];
        elseif nargin == 2
            varargout{1} = C3Ddata;
        end
        return;
    end

    % ####################
    % Adjust HSL and HSR (detected too early on swing leg retraction)
    % Shift forward until the heel velocity stops becoming more negative 
    % This doens't work when stepping on an accelerating belt
    lheeD = [0 0 0 ; diff(lhee)];
    rheeD = [0 0 0 ; diff(rhee)];
    for i = 1:length(HSL)
        n = 0;
        while ( lheeD(HSL(i,1)+n+1,2) < lheeD(HSL(i,1)+n,2) ) && ( HSL(i,1)+n+1 < size(lheeD,1) )
            n = n + 1;
        end
        HSL(i,1) = HSL(i,1) + n;
    end
    for i = 1:length(HSR)
        n = 0;
        while ( rheeD(HSR(i,1)+n+1,2) < rheeD(HSR(i,1)+n,2) ) && ( HSR(i,1)+n+1 < size(lheeD,1) )
            n = n + 1;
        end
        HSR(i,1) = HSR(i,1) + n;
    end
    % ####################
    
    % Start and end indices single supports
    SSLstart = TOR(:,1);
    SSLend = HSR(:,1)-1;
    SSRstart = TOL(:,1);
    SSRend = HSL(:,1)-1;
    
    % Remove everything that's before the first start
    if any( SSLend < SSLstart(1) )
        SSLend( 1:sum(SSLend < SSLstart(1)) ) = [];
    end
    if any( SSRend < SSRstart(1) )
        SSRend( 1:sum(SSRend < SSRstart(1)) ) = [];
    end
    
    % Check and remove for double events
    % NOTE : the second of the double event is removed, it isn't checked which is better
    doubstart = 1; doubend = 1; % Dummy assign
    while ~isempty(doubstart) && ~isempty(doubend)
        [SSLsort,idx] = sort([SSLstart ; SSLend]);
        
        indic = [ones(size(SSLstart)) ; 2.*ones(size(SSLend))];
        doubles = find( diff(indic(idx)) == 0 ) + 1;
        indic = indic(idx);
        
        doubstart = doubles(indic(doubles) == 1);
        doubend = doubles(indic(doubles) == 2);
        
        SSLstart = setdiff(SSLstart , SSLsort(doubstart));
        SSLend = setdiff(SSLend , SSLsort(doubend));
    end
    doubstart = 1; doubend = 1; % Dummy assign
    while ~isempty(doubstart) && ~isempty(doubend)
        [SSRsort,idx] = sort([SSRstart ; SSRend]);
        
        indic = [ones(size(SSRstart)) ; 2.*ones(size(SSRend))];
        doubles = find( diff(indic(idx)) == 0 ) + 1;
        indic = indic(idx);
        
        doubstart = doubles(indic(doubles) == 1);
        doubend = doubles(indic(doubles) == 2);
        
        SSRstart = setdiff(SSRstart , SSRsort(doubstart));
        SSRend = setdiff(SSRend , SSRsort(doubend));
    end

    % Remove everything that's after the last end
    if any( SSLstart > SSLend(end) )
        SSLstart( end-sum(SSLstart > SSLend(end))+1:end ) = [];
    end
    if any( SSRstart > SSRend(end) )
        SSRstart( end-sum(SSRstart > SSRend(end))+1:end ) = [];
    end
    
    % Start and end indices double support
    DSstart = sort([SSLend+1 ; SSRend+1]);
    DSend = sort([SSLstart-1 ; SSRstart-1]);

    if ~isempty(DSstart) % Not available in running
        if DSstart(end) > DSend(end)
            DSstart(end) = [];
        end
        if DSend(1) < DSstart(1)
            DSend(1) = [];
        end
        if sum(DSstart < DSend(1)) > 1  % Remove multiple starts
            DSstart(1:sum(DSstart < DSend(1))-1) = [];
        end
        if sum(DSend > DSstart(end)) > 1  % Remove multiple ends
            DSend(end - sum(DSend > DSstart(end)) +1 : end) = [];
        end
    end
    
    % Compute phases as vector
    SSL = false(size(lhee,1),1);
    for i = 1:length(SSLstart)
        SSL(SSLstart(i):SSLend(i)) = true;
    end
    SSR = false(size(lhee,1),1);
    for i = 1:length(SSRstart)
        SSR(SSRstart(i):SSRend(i)) = true;
    end
    DS = false(size(lhee,1),1);
    for i = 1:length(DSstart)
        DS(DSstart(i):DSend(i)) = true;
    end
    
    % Subphases
    % Mind that e.g. HSL does not necessarily mean the start of DSLIFOR
    % Hence we do it this way
    DSLIFOR = ( ( (lhee(:,2) - rhee(:,2)) > 0 ) + DS ) == 2;
    DSRIFOL = ( ( (rhee(:,2) - lhee(:,2)) > 0 ) + DS ) == 2;

    % Subphase onset indices
    DSLIFORstart = find(diff(DSLIFOR)>0)+1;
    DSLIFORend = find(diff(DSLIFOR)<0);
    DSRIFOLstart = find(diff(DSRIFOL)>0)+1;
    DSRIFOLend = find(diff(DSRIFOL)<0);
    
    % Remove starts without end and vice versa
    if ~isempty(DSLIFORstart)
        if DSLIFORstart(end) > DSLIFORend(end)
            DSLIFORstart(end) = [];
        end
        if DSLIFORend(1) < DSLIFORstart(1)
            DSLIFORend(1) = [];
        end
    end
    if ~isempty(DSRIFOLstart)
        if DSRIFOLstart(end) > DSRIFOLend(end)
            DSRIFOLstart(end) = [];
        end
        if DSRIFOLend(1) < DSRIFOLstart(1)
            DSRIFOLend(1) = [];
        end
    end
    
%% Store data for output
       
    phasecode = zeros(length(lhee),1);
%     phasecode(NC) = 1;  % Not available yet
    phasecode(SSL) = 2;
    phasecode(SSR) = 3;
    phasecode(DS) = 4;
    phasecode(DSLIFOR) = 5;
    phasecode(DSRIFOL) = 6;
    
    phasevector = zeros(length(lhee),7);
%     phasevector(NCend,1) = NCend - NCstart + 1;  % Not available yet
    phasevector(SSLend,2) = SSLend - SSLstart + 1;
    phasevector(SSRend,3) = SSRend - SSRstart + 1;
    phasevector(DSend,4) = DSend - DSstart + 1;
    phasevector(DSLIFORend,5) = DSLIFORend - DSLIFORstart + 1;
    phasevector(DSRIFOLend,6) = DSRIFOLend - DSRIFOLstart + 1;
    phasevector(:,7) = phasecode;
    
    phasevectorlabel = {'NC','SSL','SSR','DS','DSLIFOR','DSRIFOL','PHASECODE'};
    
   
%     phaseidx.NC = [NCstart NCend]; % Not available yet
    phaseidx.NC = [];
    phaseidx.SSL = [SSLstart SSLend];
    phaseidx.SSR = [SSRstart SSRend];
    phaseidx.DS = [DSstart DSend];
    phaseidx.DSLIFOR = [DSLIFORstart DSLIFORend];
    phaseidx.DSRIFOL = [DSRIFOLstart DSRIFOLend];
    
 
%% Generate output
    
if nargin == 1
    
    varargout{1} = phasevector;
    varargout{2} = phasevectorlabel;
    varargout{3} = phaseidx;
    
elseif nargin == 2
    
    C3Ddata.Event.GaitPhaseM = phasevector;
    C3Ddata.Event.GaitPhaseMLabel = phasevectorlabel;
    C3Ddata.Event.GaitPhaseMIdx = phaseidx;
    
    varargout{1} = C3Ddata;
end


    
end