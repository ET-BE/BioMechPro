function [varargout] = getGPf_mod(inData,fthresh)
%% Gait Phase Detection based on force data
% INPUT)
% inData : can be either:
%         - Nx12 dual force plate data OR a C3Ddata structure 
%         - C3Ddata structure
% fthresh : force threshold in N for determining contact
% 
% OUTPUT)
% If inData is a structure, output is the same structure with added fields
% 
% If inData is a matrix:
% phasevector : matrix, containing a vector for each phase.
% Non-zero indices indicate the end of each respective phase.
% The value at these indices indicates the amount of samples the phase lasted.
% These vectors are useful when calculating mean / average phase duration
% 
% Coding in phase code (last column of the phasevector):
% 1 : No Contact
% 2 : SSL
% 3 : SSR
% 4 : DS
% 5 : DSLIFOR
% 6 : DSRIFOL
% NaN : unclassified
% 
% phasevectorlabel : cell, corresponding with the columns in phasevector
% Contains the names of the phases
% 
% phaseidx : structure with the indices of the start and end of each gait phase
% 
% NOTES)
% Input is assumed an Nx12 matrix with : 
% 1-3 : Fx, Fy, Fz of the left plate
% 4-6 : Mx, My, Mz of the left plate
% 7-9 : Fx, Fy, Fz of the right plate
% 10-12 : Mx, My, Mz of the right plate
% 
% In a right handed coordinate system, where +y is the walking direction
% 
% The function assumes that the data is smooth and properly filtered !
% 
% 
% TODO
% Also make this work for an Nx6 matrix with single plate data


%% Check on input data and number of outputs

if isstruct(inData) % Input should be C3Ddata structure
    
    if nargout > 1
        error('getGPf:nout','Too many output variables specified')
    end
    
    try
        fdata = inData.Force.ForceData;
    catch
        error('getGPf:forcedata','No field Force.ForceData found in input data structure');
    end
    
else % Input is a force matrix

    if nargout > 3
        error('getGPf:nout','Too many output variables specified');
    end
    
    fdata = inData;
    
    if size(fdata,2) ~= 12
        error('getGPf:input','input data does not have 12 columns');
    end
end

% fthresh is allowed < 0, for if you somehow like to work with ground action data
    
%% Get COP

    copdata = getCOP(fdata,fthresh);
    COPL = copdata.COPL;
    COPR = copdata.COPR;

%% Get gait phases

    % Contacts
    ContactL = fdata(:,3) > fthresh; % Fz Left
    ContactR = fdata(:,9) > fthresh; % Fz Right
    if ~sum(ContactL+ContactR)
            warning(['No contact with force plate in trial ' inData.Info.Trial  '. Skipping']);
            varargout{1}=inData;
        return
    end
    % Phases
    NC = (ContactL+ContactR)==0; % No contact
    SSL = (ContactL-ContactR)==1;
    SSR = (ContactR-ContactL)==1;
    DS = (ContactL+ContactR)==2; % Note that the DS assumes that you don't stand with both legs on 1 plate.

    % Phase onset indices
    NCstart = find(diff(NC)>0)+1;
    NCend = find(diff(NC)<0);
    SSLstart = find(diff(SSL)>0)+1; % +1 to correct for the diff
    TOR=SSLstart;
    SSLend = find(diff(SSL)<0);
    HSR= SSLend;
    SSRstart = find(diff(SSR)>0)+1;
    SSRend = find(diff(SSR)<0);
    DSstart = find(diff(DS)>0)+1;
    DSend = find(diff(DS)<0);
    
    if HSR(end) > TOR(end)
        HSR(end) = [];
    end
    if TOR(1)< HSR(1) 
        TOR(1) = [];
    end
    
    
    addmanual=0;
    if addmanual
        figure
        n=[1:length(inData.Force.ForceData(:,9))];%./(inData.Force.ForceFrameRate);
        plot(n,inData.Force.ForceData(:,9),1:length(SSL),SSL*100,'*')
          %Select peaks
                prompt = {'Number'};
                title = 'Introduce missing stances';
                dims = [1 35];
                definput = {'3'};
                answer = inputdlg(prompt,title,dims,definput);
                [x,y] = ginput(answer);
                %Get sample at which starts
                i_left=str2num(x());
                %Get sample at which ends
                i_right=str2num(answer{2});
    end
    
    % Subphases
    DSLIFOR = ( ((COPL(:,2)-COPR(:,2))>0) + DS )==2; % DS Left In Front Of Right
    DSRIFOL = ( ((COPR(:,2)-COPL(:,2))>0) + DS )==2;
    
    % Subphase onset indices
    DSLIFORstart = find(diff(DSLIFOR)>0)+1;
    DSLIFORend = find(diff(DSLIFOR)<0);
    DSRIFOLstart = find(diff(DSRIFOL)>0)+1;
    DSRIFOLend = find(diff(DSRIFOL)<0);
    
    % Remove ends without start and starts without end
    if ~isempty(NCstart) && ~isempty(NCend) % Generally not available for walking
        if NCstart(end) > NCend(end)
            NCstart(end) = [];
        end
        if NCend(1) < NCstart(1)
            NCend(1) = [];
        end
    end
    
    if SSLstart(end) > SSLend(end)
        SSLstart(end) = [];
    end
    if SSLend(1) < SSLstart(1)
        SSLend(1) = [];
    end
    if SSRstart(end) > SSRend(end)
        SSRstart(end) = [];
    end
    if SSRend(1) < SSRstart(1)
        SSRend(1) = [];
    end
    
    if ~isempty(DSstart) % Not available in running
        if DSstart(end) > DSend(end)
            DSstart(end) = [];
        end
        if DSend(1) < DSstart(1)
            DSend(1) = [];
        end
    end
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
       
    phasecode = zeros(length(fdata),1);
    phasecode(NC) = 1;
    phasecode(SSL) = 2;
    phasecode(SSR) = 3;
    phasecode(DS) = 4;
    phasecode(DSLIFOR) = 5;
    phasecode(DSRIFOL) = 6;
    
    phasevector = zeros(size(fdata,1),7);
    phasevector(NCend,1) = NCend - NCstart + 1;
    phasevector(SSLend,2) = SSLend - SSLstart + 1;
    phasevector(SSRend,3) = SSRend - SSRstart + 1;
    phasevector(DSend,4) = DSend - DSstart + 1;
    phasevector(DSLIFORend,5) = DSLIFORend - DSLIFORstart + 1;
    phasevector(DSRIFOLend,6) = DSRIFOLend - DSRIFOLstart + 1;
    phasevector(:,7) = phasecode;
    
    phasevectorlabel = {'NC','SSL','SSR','DS','DSLIFOR','DSRIFOL','PHASECODE'};
    
    phaseidx.NC = [NCstart NCend];
    phaseidx.SSL = [SSLstart SSLend];
    phaseidx.SSR = [SSRstart SSRend];
    phaseidx.DS = [DSstart DSend];
    phaseidx.DSLIFOR = [DSLIFORstart DSLIFORend];
    phaseidx.DSRIFOL = [DSRIFOLstart DSRIFOLend];
    
    
%% Generate output

if isstruct(inData)
    inData.Event.GaitPhaseF = phasevector;
    inData.Event.GaitPhaseFLabel = phasevectorlabel;
    inData.Event.GaitPhaseFIdx = phaseidx;
    inData.Event.GaitRstance = [HSR TOR];
    
    varargout{1} = inData;
else
    varargout{1} = phasevector;
    varargout{2} = phasevectorlabel;
    varargout{3} = phaseidx;
end

%Plots
% figure
% plot(ContactL*700); hold on
% plot(ContactR*700)
% plot(fdata(:,3))
% plot(SSLstart,ones(1,length(SSLstart))*700,'^')
% plot(SSLend,ones(1,length(SSLend))*700,'v')
% legend({'L','R','FzR','SSLstart','SSLend'})
% 
% figure
% plot(ContactL*700); hold on
% plot(ContactR*700)
% plot(fdata(:,9))
% plot(HSR,ones(1,length(HSR))*700,'^')
% plot(TOR,ones(1,length(TOR))*700,'v')
% plot(fdata(:,3),'-')
% plot(SSRstart,ones(1,length(SSRstart))*700,'*')
% plot(SSRend,ones(1,length(SSRend))*700,'o')
% legend({'L','R','FzR','HSR','TOR','FzL','SSRstart','SSRend'})
    
end
