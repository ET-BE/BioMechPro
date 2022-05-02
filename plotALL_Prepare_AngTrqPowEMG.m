%% plotAll_MF_MeanBased
% 
% Process all data, finding means

% ##################################
% SET 1 : ML SLOW
% SET 2 : ML FAST
% SET 3 : AP SLOW
% SET 4 : AP FAST
% ##################################

clear variables; close all; clc;

%% Set info

subjects = 1:10;  % Dropped subject 3 in piece collection

g = 9.81;

col = jet(8); col(1:4,:) = col(4:-1:1,:); % Swap colors of first four
colb = [0.75 0.75 0.75];
lin1 = {'<','<','<','<','>','>','>','>'};
lin2 = {'v','v','v','v','^','^','^','^'};
linw = 0.5;
fnz = 10;
mrksiz1 = 3;

%% Pre alloc 
% (only partially, as repetitions may vary per subject)

nsmpl = 50;

% Original time series
eventDatab = zeros(300,1,4,10); % data irep iset isubj (zeros, not NaN!!)
jAngDatab = NaN(300,9,3,1,4,10);
jAngDataDb = NaN(300,9,3,1,4,10);
jTrqDatab = NaN(300,9,3,1,4,10);
% comDataDb = NaN(300,10,3,1,4,10);
emgDatab = NaN(3000,12,1,4,10);
% fDatab = NaN(3000,12,1,4,10);

eventDatap = zeros(300,1,8,4,10); % data irep ipert iset isubj
jAngDatap = NaN(300,9,3,1,8,4,10);
jAngDataDp = NaN(300,9,3,1,8,4,10);
jTrqDatap = NaN(300,9,3,1,8,4,10);
% comDataDp = NaN(300,10,3,1,8,4,10);
emgDatap = NaN(3000,12,1,8,4,10);
% fDatap = NaN(3000,12,1,8,4,10);

% Data in resamples sequence
% NOTE : THE LAST NUMBER IS THE NUMBER OF SEQUENCES (5): TOR-HSR, HSR-TOL, TOL-HSL, HSL-TOR, TOR-TOR
nseq = 7;
jAngDatabt = NaN(nsmpl,9,3,1,4,10,nseq);
jAngDataDbt = NaN(nsmpl,9,3,1,4,10,nseq);
jTrqDatabt = NaN(nsmpl,9,3,1,4,10,nseq);
jPowDatabt = NaN(nsmpl,9,3,1,4,10,nseq);
% jEgyNegDatabt = NaN(nsmpl,9,3,1,4,10,nseq);
% jEgyPosDatabt = NaN(nsmpl,9,3,1,4,10,nseq);
% comDataDbt = NaN(nsmpl,10,3,1,4,10,nseq);
emgDatabt = NaN(nsmpl,12,1,4,10,nseq);
% fDatabt = NaN(nsmpl,12,1,4,10,nseq);

jAngDatapt = NaN(nsmpl,9,3,1,8,4,10,nseq);
jAngDataDpt = NaN(nsmpl,9,3,1,8,4,10,nseq);
jTrqDatapt = NaN(nsmpl,9,3,1,8,4,10,nseq);
jPowDatapt = NaN(nsmpl,9,3,1,8,4,10,nseq);
% jEgyNegDatapt = NaN(nsmpl,9,3,1,8,4,10,nseq);
% jEgyPosDatapt = NaN(nsmpl,9,3,1,8,4,10,nseq);
% comDataDpt = NaN(nsmpl,10,3,1,8,4,10,nseq);
emgDatapt = NaN(nsmpl,12,1,8,4,10,nseq);
% fDatapt = NaN(nsmpl,12,1,8,4,10,nseq);

% Data in events
nevt = 6;
tDatabe = NaN(nevt,1,4,10);
tDatape = NaN(nevt,1,8,4,10);

% Scaled constants
% w0 = NaN(4,10); % XCOM variable

%% Collecting data

% Collect data
for isubj = subjects
    
    getInfo;
    
    % Load Marker piece
    load([subjroot '\pieceMarker.mat']);
    
    % Load Event piece
    load([subjroot '\pieceEvent.mat']);
    
    % Load EMG piece
    load([subjroot '\pieceEMG.mat']);
    
    % Load force piece (Tycho)
%     load([subjroot '\pieceForce.mat']);
    
    % Sample freq
    fs_mrk = pieceMarker.Set(1).Pert(1).Data.MarkerFrameRate;
    fs_emg = pieceEMG.Set(1).Pert(1).Data.EMGFrameRate;

    % Loop through the sets
    missed = 0;
    for iset = 1:4
        
        % #####################################
        % ##            BASELINE             ##
        % #####################################
        % Time series
        eHSR = squeeze(pieceEvent.Set(iset).Base(1).Data.GaitPhaseM(:,2,:)); 
        eTOL = squeeze(pieceEvent.Set(iset).Base(1).Data.GaitPhaseM(:,6,:)); 
        eHSL = squeeze(pieceEvent.Set(iset).Base(1).Data.GaitPhaseM(:,3,:)); 
        eTOR = squeeze(pieceEvent.Set(iset).Base(1).Data.GaitPhaseM(:,5,:)); 
        event = eHSR + eTOL + eHSL + eTOR;
        event(1:65,:) = 0; % Remove all before TOR, so we know order is correct (start with HSR)
        event([50 65],:) = 1; % Insert start of sequence at TOR (and fictional end)
        eventDatab(:,1:size(event,2),iset,isubj) = event;
        
        dat = pieceMarker.Set(iset).Base(1).Data.JointAngData;
        jAngDatab(:,:,:,1:size(dat,4),iset,isubj) = dat;

        dat = pieceMarker.Set(iset).Base(1).Data.JointAngDataD;
        jAngDataDb(:,:,:,1:size(dat,4),iset,isubj) = dat;
        
        dat = pieceMarker.Set(iset).Base(1).Data.JointTrqData;
        jTrqDatab(:,:,:,1:size(dat,4),iset,isubj) = dat;
        
%         dat = pieceMarker.Set(iset).Base(1).Data.COMDataD;
%         comDataDb(:,:,:,1:size(dat,4),iset,isubj) = dat;
        
        dat = pieceEMG.Set(iset).Base(1).Data.EMGData;
        emgDatab(:,:,1:size(dat,3),iset,isubj) = dat;
        
%         dat = pieceForce.Set(iset).Base(1).Data.ForceData;
%         fDatab(:,:,1:size(dat,3),iset,isubj) = dat;
        
        % #####################################
        % ##          PERTURBATIONS          ##
        % #####################################
        for ipert = 1:8
            eHSR = squeeze(pieceEvent.Set(iset).Pert(ipert).Data.GaitPhaseM(:,2,:)); 
            eTOL = squeeze(pieceEvent.Set(iset).Pert(ipert).Data.GaitPhaseM(:,6,:)); 
            eHSL = squeeze(pieceEvent.Set(iset).Pert(ipert).Data.GaitPhaseM(:,3,:)); 
            eTOR = squeeze(pieceEvent.Set(iset).Pert(ipert).Data.GaitPhaseM(:,5,:)); 
            event = eHSR + eTOL + eHSL + eTOR;
            event(1:65,:) = 0; % Remove all before TOR, so we know order is correct (start with HSR)
            event([50 65],:) = 1; % Insert start and end of perturbation
            eventDatap(:,1:size(event,2),ipert,iset,isubj) = event;
            
            dat = pieceMarker.Set(iset).Pert(ipert).Data.JointAngData;
            jAngDatap(:,:,:,1:size(dat,4),ipert,iset,isubj) = dat;
            
            dat = pieceMarker.Set(iset).Pert(ipert).Data.JointAngDataD;
            jAngDataDp(:,:,:,1:size(dat,4),ipert,iset,isubj) = dat;
            
            dat = pieceMarker.Set(iset).Pert(ipert).Data.JointTrqData;
            jTrqDatap(:,:,:,1:size(dat,4),ipert,iset,isubj) = dat;
            
%             dat = pieceMarker.Set(iset).Pert(ipert).Data.COMDataD;
%             comDataDp(:,:,:,1:size(dat,4),ipert,iset,isubj) = dat;
            
            dat = pieceEMG.Set(iset).Pert(ipert).Data.EMGData;
            emgDatap(:,:,1:size(dat,3),ipert,iset,isubj) = dat;
            
%             dat = pieceForce.Set(iset).Pert(ipert).Data.ForceData;
%             fDatap(:,:,1:size(dat,3),ipert,iset,isubj) = dat;
            
        end
    end
end

%% Free some memory
clear pieceMarker pieceForce pieceEvent foo dat
clear pieceForce

%% Make filler zeros equal to NaN

jAngDatab(jAngDatab == 0) = NaN;
jAngDataDb(jAngDataDb == 0) = NaN;
jTrqDatab(jTrqDatab == 0) = NaN;
% comDataDb(comDataDb == 0) = NaN;
emgDatab(emgDatab == 0) = NaN;
% fDatab(fDatab == 0) = NaN;

jAngDatap(jAngDatap == 0) = NaN;
jAngDataDp(jAngDataDp == 0) = NaN;
jTrqDatap(jTrqDatap == 0) = NaN;
% comDataDp(comDataDp == 0) = NaN;
emgDatap(emgDatap == 0) = NaN;
% fDatap(fDatap == 0) = NaN;

%% Downsample EMG data

% Not sure yet if I want to do this.
% But required if you use events based on markerdata
emgDatab = downsample(emgDatab,10);
emgDatap = downsample(emgDatap,10);

% fDatab = downsample(fDatab,10);
% fDatap = downsample(fDatap,10);

%% Get joint power

jPowDatab = jAngDataDb.*jTrqDatab;
jPowDatap = jAngDataDp.*jTrqDatap;

%% Remove outliers
%
% NOTE: THERE ARE SOME MAJOR ISSUES WITH SUBJECT 10 ANGLES & TORQUES IN SET 2 
% (fast, ML)

% Joint Angles & velocities
prct = 90;
outSmpl = 100;
outMaskA = sum( abs(jAngDatab - repmat(nanmedian(jAngDatab,4) ,  [1 1 1 size(jAngDatab,4) 1 1])) > repmat( prctile( abs(jAngDatab - repmat(nanmedian(jAngDatab,4) , [1 1 1 size(jAngDatab,4) 1 1])) , prct , 4 ), [1 1 1 size(jAngDatab,4) 1 1]) ) > outSmpl;
outMaskB = sum(abs(jAngDatab)>0.5*pi,1) > 5;
outMaskC = sum( abs(jAngDataDb - repmat(nanmedian(jAngDataDb,4) ,  [1 1 1 size(jAngDataDb,4) 1 1])) > repmat( prctile( abs(jAngDataDb - repmat(nanmedian(jAngDataDb,4) , [1 1 1 size(jAngDataDb,4) 1 1])) , prct , 4 ), [1 1 1 size(jAngDataDb,4) 1 1]) ) > outSmpl;
outMaskD = sum(abs(jAngDataDb)>11,1) > 1;
outMask = outMaskA | outMaskB | outMaskC | outMaskD;
outMask = repmat(outMask, [size(jAngDatab,1) 1 1 1 1 1]);
% foo = jAngDataDb;
% foo(outMask) = NaN;
% return;
jAngDatab(outMask) = NaN;
jAngDataDb(outMask) = NaN;

prct = 99; % Repeat once more for large but short duration deviations
outSmpl = 10;
outMaskA = sum( abs(jAngDatab - repmat(nanmedian(jAngDatab,4) ,  [1 1 1 size(jAngDatab,4) 1 1])) > repmat( prctile( abs(jAngDatab - repmat(nanmedian(jAngDatab,4) , [1 1 1 size(jAngDatab,4) 1 1])) , prct , 4 ), [1 1 1 size(jAngDatab,4) 1 1]) ) > outSmpl;
outMaskC = sum( abs(jAngDataDb - repmat(nanmedian(jAngDataDb,4) ,  [1 1 1 size(jAngDataDb,4) 1 1])) > repmat( prctile( abs(jAngDataDb - repmat(nanmedian(jAngDataDb,4) , [1 1 1 size(jAngDataDb,4) 1 1])) , prct , 4 ), [1 1 1 size(jAngDataDb,4) 1 1]) ) > outSmpl;
outMask = outMaskA | outMaskC;
outMask = repmat(outMask, [size(jAngDatab,1) 1 1 1 1 1]);
jAngDatab(outMask) = NaN;
jAngDataDb(outMask) = NaN;

% Remove subject 10 in set 2
jAngDatab(:,:,:,:,2,10) = NaN;
jAngDataDb(:,:,:,:,2,10) = NaN;

prct = 90;
outSmpl = 100;
outMaskA = sum( abs(jAngDatap - repmat(nanmedian(jAngDatap,4) ,  [1 1 1 size(jAngDatap,4) 1 1 1])) > repmat( prctile( abs(jAngDatap - repmat(nanmedian(jAngDatap,4) , [1 1 1 size(jAngDatap,4) 1 1 1])) , prct , 4 ), [1 1 1 size(jAngDatap,4) 1 1 1]) ) > outSmpl;
outMaskB = sum(abs(jAngDatap)> 0.75*pi,1) > 5;
outMaskC = sum( abs(jAngDataDp - repmat(nanmedian(jAngDataDp,4) ,  [1 1 1 size(jAngDataDp,4) 1 1 1])) > repmat( prctile( abs(jAngDataDp - repmat(nanmedian(jAngDataDp,4) , [1 1 1 size(jAngDataDp,4) 1 1 1])) , prct , 4 ), [1 1 1 size(jAngDataDp,4) 1 1 1]) ) > outSmpl;
outMaskD = sum(abs(jAngDataDp)> 11,1) > 1;
outMask = outMaskA | outMaskB | outMaskC | outMaskD;
outMask = repmat(outMask, [size(jAngDatap,1) 1 1 1 1 1]);
% foo = jAngDatap;
% bar = jAngDataDp;
% foo(outMask) = NaN;
% bar(outMask) = NaN;
% return;
jAngDatap(outMask) = NaN;
jAngDataDp(outMask) = NaN;

% prct = 99;  % Repeat once more for large but short duration deviations
% outSmpl = 10;
% outMaskA = sum( abs(jAngDatap - repmat(nanmedian(jAngDatap,4) ,  [1 1 1 size(jAngDatap,4) 1 1 1])) > repmat( prctile( abs(jAngDatap - repmat(nanmedian(jAngDatap,4) , [1 1 1 size(jAngDatap,4) 1 1 1])) , prct , 4 ), [1 1 1 size(jAngDatap,4) 1 1 1]) ) > outSmpl;
% outMaskC = sum( abs(jAngDataDp - repmat(nanmedian(jAngDataDp,4) ,  [1 1 1 size(jAngDataDp,4) 1 1 1])) > repmat( prctile( abs(jAngDataDp - repmat(nanmedian(jAngDataDp,4) , [1 1 1 size(jAngDataDp,4) 1 1 1])) , prct , 4 ), [1 1 1 size(jAngDataDp,4) 1 1 1]) ) > outSmpl;
% outMask = outMaskA | outMaskC;
% outMask = repmat(outMask, [size(jAngDatap,1) 1 1 1 1 1]);
% jAngDatap(outMask) = NaN;
% jAngDataDp(outMask) = NaN;

% Remove subject 10 in set 2
jAngDatap(:,:,:,:,:,2,10) = NaN;
jAngDataDp(:,:,:,:,:,2,10) = NaN;

% Joint Torques
prct = 90;
outSmpl = 100;
outMaskA = sum( abs(jTrqDatab - repmat(nanmedian(jTrqDatab,4) ,  [1 1 1 size(jTrqDatab,4) 1 1])) > repmat( prctile( abs(jTrqDatab - repmat(nanmedian(jTrqDatab,4) , [1 1 1 size(jTrqDatab,4) 1 1])) , prct , 4 ), [1 1 1 size(jTrqDatab,4) 1 1]) ) > outSmpl;
outMaskB = sum(abs(jTrqDatab)>170,1) > 5;
outMask = outMaskA | outMaskB;
outMask = repmat(outMask, [size(jTrqDatab,1) 1 1 1 1 1]);
% foo = jTrqDatab;
% foo(outMask) = NaN;
% return;
jTrqDatab(outMask) = NaN;

prct = 99; % Repeat once more for large but short duration deviations
outSmpl = 10;
outMaskA = sum( abs(jTrqDatab - repmat(nanmedian(jTrqDatab,4) ,  [1 1 1 size(jTrqDatab,4) 1 1])) > repmat( prctile( abs(jTrqDatab - repmat(nanmedian(jTrqDatab,4) , [1 1 1 size(jTrqDatab,4) 1 1])) , prct , 4 ), [1 1 1 size(jTrqDatab,4) 1 1]) ) > outSmpl;
outMask = repmat(outMaskA, [size(jTrqDatab,1) 1 1 1 1 1]);
jTrqDatab(outMask) = NaN;

% Remove subject 10 in set 2
jTrqDatab(:,:,:,:,2,10) = NaN;

prct = 90;
outSmpl = 100;
outMaskA = sum( abs(jTrqDatap - repmat(nanmedian(jTrqDatap,4) ,  [1 1 1 size(jTrqDatap,4) 1 1])) > repmat( prctile( abs(jTrqDatap - repmat(nanmedian(jTrqDatap,4) , [1 1 1 size(jTrqDatap,4) 1 1])) , prct , 4 ), [1 1 1 size(jTrqDatap,4) 1 1]) ) > outSmpl;
outMaskB = sum(abs(jTrqDatap)>170,1) > 5;
outMask = outMaskA | outMaskB;
outMask = repmat(outMask, [size(jTrqDatap,1) 1 1 1 1 1]);
% foo = jTrqDatap;
% foo(outMask) = NaN;
jTrqDatap(outMask) = NaN;
% return;
% prct = 99; % Repeat once more for large but short duration deviations
% outSmpl = 10;
% outMaskA = sum( abs(jTrqDatap - repmat(nanmedian(jTrqDatap,4) ,  [1 1 1 size(jTrqDatap,4) 1 1])) > repmat( prctile( abs(jTrqDatap - repmat(nanmedian(jTrqDatap,4) , [1 1 1 size(jTrqDatap,4) 1 1])) , prct , 4 ), [1 1 1 size(jTrqDatap,4) 1 1]) ) > outSmpl;
% outMask = repmat(outMaskA, [size(jTrqDatap,1) 1 1 1 1 1]);
% foo = jTrqDatap;
% foo(outMask) = NaN;
% % jTrqDatap(outMask) = NaN;

% Remove subject 10 in set 2
jTrqDatap(:,:,:,:,2,10) = NaN;


% Joint Power
prct = 90;
outSmpl = 100;
outMaskA = sum( abs(jPowDatab - repmat(nanmedian(jPowDatab,4) ,  [1 1 1 size(jPowDatab,4) 1 1])) > repmat( prctile( abs(jPowDatab - repmat(nanmedian(jPowDatab,4) , [1 1 1 size(jPowDatab,4) 1 1])) , prct , 4 ), [1 1 1 size(jPowDatab,4) 1 1]) ) > outSmpl;
outMaskB = sum(abs(jPowDatab)>170,1) > 5;
outMask = outMaskA | outMaskB;
outMask = repmat(outMask, [size(jPowDatab,1) 1 1 1 1 1]);
% foo = jPowDatab;
% foo(outMask) = NaN;
jPowDatab(outMask) = NaN;

prct = 99; % Repeat once more for large but short duration deviations
outSmpl = 10;
outMaskA = sum( abs(jPowDatab - repmat(nanmedian(jPowDatab,4) ,  [1 1 1 size(jPowDatab,4) 1 1])) > repmat( prctile( abs(jPowDatab - repmat(nanmedian(jPowDatab,4) , [1 1 1 size(jPowDatab,4) 1 1])) , prct , 4 ), [1 1 1 size(jPowDatab,4) 1 1]) ) > outSmpl;
outMask = repmat(outMaskA, [size(jPowDatab,1) 1 1 1 1 1]);
jPowDatab(outMask) = NaN;

prct = 90;
outSmpl = 100;
outMaskA = sum( abs(jPowDatap - repmat(nanmedian(jPowDatap,4) ,  [1 1 1 size(jPowDatap,4) 1 1])) > repmat( prctile( abs(jPowDatap - repmat(nanmedian(jPowDatap,4) , [1 1 1 size(jPowDatap,4) 1 1])) , prct , 4 ), [1 1 1 size(jPowDatap,4) 1 1]) ) > outSmpl;
outMaskB = sum(abs(jPowDatap)>170,1) > 5;
outMask = outMaskA | outMaskB;
outMask = repmat(outMask, [size(jPowDatap,1) 1 1 1 1 1]);
% foo = jPowDatap;
% foo(outMask) = NaN;
jPowDatap(outMask) = NaN;


% EMG
prct = 90;
outSmpl = 100;
outMaskA = sum( abs(emgDatab - repmat(nanmedian(emgDatab,3) ,  [1 1 size(emgDatab,3) 1 1])) > repmat( prctile( abs(emgDatab - repmat(nanmedian(emgDatab,3) , [1 1 size(emgDatab,3) 1 1])) , prct , 3 ), [1 1 size(emgDatab,3) 1 1]) ) > outSmpl;
outMask = outMaskA;
outMask = repmat(outMask, [size(emgDatab,1) 1 1 1 1 1]);
% foo = emgDatab;
% foo(outMask) = NaN;
% return;
emgDatab(outMask) = NaN;

prct = 90;
outSmpl = 100;
outMaskA = sum( abs(emgDatap - repmat(nanmedian(emgDatap,3) ,  [1 1 size(emgDatap,3) 1 1 1])) > repmat( prctile( abs(emgDatap - repmat(nanmedian(emgDatap,3) , [1 1 size(emgDatap,3) 1 1 1])) , prct , 3 ), [1 1 size(emgDatap,3) 1 1 1]) ) > outSmpl;
outMask = outMaskA;
outMask = repmat(outMask, [size(emgDatap,1) 1 1 1 1 1]);
% foo = emgDatab;
% foo(outMask) = NaN;
% return;
emgDatap(outMask) = NaN;

% Remove subject 10 adductor longi (Clipping Adductor Longi)
eDatab(:,[6 12],:,:,10) = NaN;
eDatap(:,[6 12],:,:,:,10) = NaN;


%% Get integrated joint power

% foo = jPowDatab; foo(foo>0) = 0;
% bar = jPowDatab; bar(bar<0) = 0;
% jEgyNegDatab = cumtrapz( foo , 1) ./ fs_mrk;
% jEgyPosDatab = cumtrapz( bar , 1) ./ fs_mrk;
% 
% foo = jPowDatap; foo(foo>0) = 0;
% bar = jPowDatap; bar(bar<0) = 0;
% jEgyNegDatap = cumtrapz( foo , 1) ./ fs_mrk;
% jEgyPosDatap = cumtrapz( bar , 1) ./ fs_mrk;


%% Add belt velocity to global COM velocity

% for isubj = subjects
%     
%     getInfo;
% 
%     comDataDb(:,end,2,:,[1 3],isubj) = comDataDb(:,end,2,:,[1 3],isubj) + (2.25/3.6)*sqrt(subjleglen);
%     comDataDp(:,end,2,:,:,[1 3],isubj) = comDataDp(:,end,2,:,:,[1 3],isubj) + (2.25/3.6).*sqrt(subjleglen);
% 
%     comDataDb(:,end,2,:,[2 4],isubj) = comDataDb(:,end,2,:,[2 4],isubj) + (4.5/3.6)*sqrt(subjleglen);
%     comDataDp(:,end,2,:,:,[2 4],isubj) = comDataDp(:,end,2,:,:,[2 4],isubj) + (4.5/3.6).*sqrt(subjleglen);
% 
% end

%% Get resampled time series

for isubj = subjects
    for iset = 1:4

        % Event based
        foo = getEventData(eventDatab(:,:,iset,isubj),eventDatab(:,:,iset,isubj),1,2,nevt) ./ fs_mrk;
        tDatabe(:,1:size(foo,2),iset,isubj) = foo;
        
        % Time sequence based
        for iseq = 1:6
            switch iseq
                case 1
                    evts = [1 2]; % TOR-STOP
                case 2
                    evts = [2 3]; % STOP-HSR
                case 3
                    evts = [1 3]; % TOR-HSR
                case 4
                    evts = [3 4]; % HSR-TOL
                case 5
                    evts = [4 5]; % TOL-HSL
                case 6
                    evts = [5 6]; % HSL-TOR
                case 7
                    evts = [1 6]; % TOR-TOR
            end
            foo = getEventSequence(jAngDatab(:,:,:,:,iset,isubj),eventDatab(:,:,iset,isubj),1,4,nsmpl,evts);
            jAngDatabt(:,:,:,1:size(foo,4),iset,isubj,iseq) = foo;

            foo = getEventSequence(jAngDataDb(:,:,:,:,iset,isubj),eventDatab(:,:,iset,isubj),1,4,nsmpl,evts);
            jAngDataDbt(:,:,:,1:size(foo,4),iset,isubj,iseq) = foo;
            
            foo = getEventSequence(jTrqDatab(:,:,:,:,iset,isubj),eventDatab(:,:,iset,isubj),1,4,nsmpl,evts);
            jTrqDatabt(:,:,:,1:size(foo,4),iset,isubj,iseq) = foo;
            
            foo = getEventSequence(jPowDatab(:,:,:,:,iset,isubj),eventDatab(:,:,iset,isubj),1,4,nsmpl,evts);
            jPowDatabt(:,:,:,1:size(foo,4),iset,isubj,iseq) = foo;
            
%             foo = getEventSequence(jEgyNegDatab(:,:,:,:,iset,isubj),eventDatab(:,:,iset,isubj),1,4,nsmpl,evts);
%             jEgyNegDatabt(:,:,:,1:size(foo,4),iset,isubj,iseq) = foo;
%             
%             foo = getEventSequence(jEgyPosDatab(:,:,:,:,iset,isubj),eventDatab(:,:,iset,isubj),1,4,nsmpl,evts);
%             jEgyPosDatabt(:,:,:,1:size(foo,4),iset,isubj,iseq) = foo;
%             
%             foo = getEventSequence(comDataDb(:,:,:,:,iset,isubj),eventDatab(:,:,iset,isubj),1,4,nsmpl,evts);
%             comDataDbt(:,:,:,1:size(foo,4),iset,isubj,iseq) = foo;
            
            foo = getEventSequence(emgDatab(:,:,:,iset,isubj),eventDatab(:,:,iset,isubj),1,3,nsmpl,evts);
            emgDatabt(:,:,1:size(foo,3),iset,isubj,iseq) = foo;
            
%             foo = getEventSequence(fDatab(:,:,:,iset,isubj),eventDatab(:,:,iset,isubj),1,3,nsmpl,evts);
%             fDatabt(:,:,1:size(foo,3),iset,isubj,iseq) = foo;
            
        end
        
        for ipert = 1:8
            
            % Event based
            foo = getEventData(eventDatap(:,:,ipert,iset,isubj),eventDatap(:,:,ipert,iset,isubj),1,2,nevt) ./ fs_mrk;
            tDatape(:,1:size(foo,2),ipert,iset,isubj) = foo;
            
            % Time sequence based
            for iseq = 1:6
                switch iseq
                    case 1
                        evts = [1 2]; % Start(TOR)-STOP
                    case 2
                        evts = [2 3]; % STOP-HSR
                    case 3
                        evts = [1 3]; % Start-HSR
                    case 4
                        evts = [3 4]; % HSR-TOL
                    case 5
                        evts = [4 5]; % TOL-HSL
                    case 6
                        evts = [5 6]; % HSL-TOR
                    case 7
                        evts = [1 6]; % TOR-TOR
                end
                
                foo = getEventSequence(jAngDatap(:,:,:,:,ipert,iset,isubj),eventDatap(:,:,ipert,iset,isubj),1,4,nsmpl,evts);
                jAngDatapt(:,:,:,1:size(foo,4),ipert,iset,isubj,iseq) = foo;
                
                foo = getEventSequence(jAngDataDp(:,:,:,:,ipert,iset,isubj),eventDatap(:,:,ipert,iset,isubj),1,4,nsmpl,evts);
                jAngDataDpt(:,:,:,1:size(foo,4),ipert,iset,isubj,iseq) = foo;
                
                foo = getEventSequence(jTrqDatap(:,:,:,:,ipert,iset,isubj),eventDatap(:,:,ipert,iset,isubj),1,4,nsmpl,evts);
                jTrqDatapt(:,:,:,1:size(foo,4),ipert,iset,isubj,iseq) = foo;
                
                foo = getEventSequence(jPowDatap(:,:,:,:,ipert,iset,isubj),eventDatap(:,:,ipert,iset,isubj),1,4,nsmpl,evts);
                jPowDatapt(:,:,:,1:size(foo,4),ipert,iset,isubj,iseq) = foo;

%                 foo = getEventSequence(jEgyNegDatap(:,:,:,:,ipert,iset,isubj),eventDatap(:,:,ipert,iset,isubj),1,4,nsmpl,evts);
%                 jEgyNegDatabt(:,:,:,1:size(foo,4),ipert,iset,isubj,iseq) = foo;
% 
%                 foo = getEventSequence(jEgyPosDatap(:,:,:,:,ipert,iset,isubj),eventDatap(:,:,ipert,iset,isubj),1,4,nsmpl,evts);
%                 jEgyPosDatabt(:,:,:,1:size(foo,4),ipert,iset,isubj,iseq) = foo;
%                 
%                 foo = getEventSequence(comDataDp(:,:,:,:,ipert,iset,isubj),eventDatap(:,:,ipert,iset,isubj),1,4,nsmpl,evts);
%                 comDataDpt(:,:,:,1:size(foo,4),ipert,iset,isubj,iseq) = foo;
%                 
                foo = getEventSequence(emgDatap(:,:,:,ipert,iset,isubj),eventDatap(:,:,ipert,iset,isubj),1,3,nsmpl,evts);
                emgDatapt(:,:,1:size(foo,3),ipert,iset,isubj,iseq) = foo;
                
%                 foo = getEventSequence(fDatap(:,:,:,ipert,iset,isubj),eventDatap(:,:,ipert,iset,isubj),1,3,nsmpl,evts);
%                 fDatapt(:,:,1:size(foo,3),ipert,iset,iset,iseq) = foo;
                
            end
        end
    end
end

% Clean up original storage matrices
clear eventDatab jAngDatab jAngDataDb jTrqDatab jPowDatab jEgyNegDatab jEgyPosDatab comDataDb emgDatab ...
    eventDatap jAngDatap jAngDataDp jTrqDatap jPowDatap jEgyNegDatap jEgyPosDatap comDataDp emgDatap;
% clear fDatab fDatap


%% Make filler zeros equal to NaN, again
% Zero data can occur due "filling" of matrices
% (except for in eventData, as here the zero's are required

jAngDatabt(jAngDatabt == 0) = NaN;
jAngDataDbt(jAngDataDbt == 0) = NaN;
jTrqDatabt(jTrqDatabt == 0) = NaN;
jPowDatabt(jPowDatabt == 0) = NaN;
% comDataDbt(comDataDbt == 0) = NaN;
emgDatabt(emgDatabt == 0) = NaN;
% fDatabt(fDatabt == 0) = NaN;

jAngDatapt(jAngDatapt == 0) = NaN;
jAngDataDpt(jAngDataDpt == 0) = NaN;
jTrqDatapt(jTrqDatapt == 0) = NaN;
jPowDatapt(jPowDatapt == 0) = NaN;
% comDataDpt(comDatDpt == 0) = NaN;
emgDatapt(emgDatapt == 0) = NaN;
% fDatapt(fDatapt == 0) = NaN;

tDatabe(tDatabe == 0) = NaN;
tDatape(tDatape == 0) = NaN;

%% Do scaling
% 
% Note: you don't have to scale angles, as these are already dimensionless

% Do scaling
for isubj = subjects
    
    getInfo;
    
    for iset = 1:4

        % Scale to body weight and height
        jTrqDatabt(:,:,:,:,iset,isubj,:) = jTrqDatabt(:,:,:,:,iset,isubj,:) ./ (subjmass*g*subjheight);
        jTrqDatapt(:,:,:,:,:,iset,isubj,:) = jTrqDatapt(:,:,:,:,:,iset,isubj,:) ./ (subjmass*g*subjheight);
        
        jPowDatabt(:,:,:,:,iset,isubj,:) = jPowDatabt(:,:,:,:,iset,isubj,:) ./ (subjmass*sqrt(g)*(subjheight.^1.5));
        jPowDatapt(:,:,:,:,:,iset,isubj,:) = jPowDatapt(:,:,:,:,:,iset,isubj,:) ./ (subjmass*sqrt(g)*(subjheight.^1.5));
        
%         jEgyNegDatabt(:,:,:,:,iset,isubj,:) = jEgyNegDatabt(:,:,:,:,iset,isubj,:) ./ (subjmass*g*subjheight);
%         jEgyPosDatabt(:,:,:,:,iset,isubj,:) = jEgyPosDatabt(:,:,:,:,iset,isubj,:) ./ (subjmass*g*subjheight);
% 
%         jEgyNegDatapt(:,:,:,:,:,iset,isubj,:) = jEgyNegDatapt(:,:,:,:,:,iset,isubj,:) ./ (subjmass*g*subjheight);
%         jEgyPosDatapt(:,:,:,:,:,iset,isubj,:) = jEgyPosDatapt(:,:,:,:,:,iset,isubj,:) ./ (subjmass*g*subjheight);
                

        % COM velocity (requires COM position and feet positions)
        
    end

%     fDatabt(:,:,:,:,isubj,:) = fDatabt(:,:,:,:,isubj,:) / subjmass;
%     fDatapt(:,:,:,:,:,isubj,:) = fDatapt(:,:,:,:,:,isubj,:) / subjmass;
    
end

% EMG scaling to individual muscles, median of max values in the baseline.
% Scale all sequences of the same muscle the same
scalemask = nanmedian(max(max(emgDatabt,[],1),[],6),3);
emgDatabt = emgDatabt ./ repmat(scalemask,[size(emgDatabt,1) 1 size(emgDatabt,3) 1 1 nseq]);
emgDatapt = emgDatapt ./ permute( repmat(scalemask,[size(emgDatapt,1) 1 size(emgDatapt,3) 1 1 nseq 8]) , [1:3 7 4:6]);

%% Get repetition averages

% Averages
jAngDatabt_m1 = squeeze(nanmean(jAngDatabt,4));
jAngDataDbt_m1 = squeeze(nanmean(jAngDataDbt,4));
jTrqDatabt_m1 = squeeze(nanmean(jTrqDatabt,4));
jPowDatabt_m1 = squeeze(nanmean(jPowDatabt,4));
% jEgyNegDatabt_m1 = squeeze(nanmean(jEgyNegDatabt,4));
% jEgyPosDatabt_m1 = squeeze(nanmean(jEgyPosDatabt,4));
% comDataDbt_m1 = squeeze(nanmean(comDataDbt,4));
emgDatabt_m1 = squeeze(nanmean(emgDatabt,3));
% fDatabt_m1 = squeeze(nanmean(fDatabt,3));

jAngDatapt_m1 = squeeze(nanmean(jAngDatapt,4));
jAngDataDpt_m1 = squeeze(nanmean(jAngDataDpt,4));
jTrqDatapt_m1 = squeeze(nanmean(jTrqDatapt,4));
jPowDatapt_m1 = squeeze(nanmean(jPowDatapt,4));
% jEgyNegDatapt_m1 = squeeze(nanmean(jEgyNegDatapt,4));
% jEgyPosDatapt_m1 = squeeze(nanmean(jEgyPosDatapt,4));
% comDataDpt_m1 = squeeze(nanmean(comDataDpt,4));
emgDatapt_m1 = squeeze(nanmean(emgDatapt,3));
% fDatapt_m1 = squeeze(nanmean(fDatapt,3));

tDatabe_m1 = squeeze(nanmean(tDatabe,2));
tDatape_m1 = squeeze(nanmean(tDatape,2));

% Standard deviations
jAngDatabt_std1 = squeeze(nanstd(jAngDatabt,[],4));
jAngDataDbt_std1 = squeeze(nanstd(jAngDataDbt,[],4));
jTrqDatabt_std1 = squeeze(nanstd(jTrqDatabt,[],4));
jPowDatabt_std1 = squeeze(nanstd(jPowDatabt,[],4));
% comDataDbt_std1 = squeeze(nanstd(comDataDbt,[],4));
emgDatabt_std1 = squeeze(nanstd(emgDatabt,[],3));
% 
jAngDatapt_std1 = squeeze(nanstd(jAngDatapt,[],4));
jAngDataDpt_std1 = squeeze(nanstd(jAngDataDpt,[],4));
jTrqDatapt_std1 = squeeze(nanstd(jTrqDatapt,[],4));
jPowDatapt_std1 = squeeze(nanstd(jPowDatapt,[],4));
% comDataDpt_std1 = squeeze(nanstd(comDataDpt,[],4));
emgDatapt_std1 = squeeze(nanstd(emgDatapt,[],3));

%% Get subject averages

% Averages
jAngDatabt_m2 = squeeze(nanmean(jAngDatabt_m1,5));
jAngDataDbt_m2 = squeeze(nanmean(jAngDataDbt_m1,5));
jTrqDatabt_m2 = squeeze(nanmean(jTrqDatabt_m1,5));
jPowDatabt_m2 = squeeze(nanmean(jPowDatabt_m1,5));
% jEgyNegDatabt_m2 = squeeze(nanmean(jEgyNegDatabt_m1,5));
% jEgyPosDatabt_m2 = squeeze(nanmean(jEgyPosDatabt_m1,5));
% comDataDbt_m2 = squeeze(nanmean(comDataDbt_m1,5));
emgDatabt_m2 = squeeze(nanmean(emgDatabt_m1,4));
% fDatabt_m2 = squeeze(nanmean(fDatabt_m1,4));

jAngDatapt_m2 = squeeze(nanmean(jAngDatapt_m1,6));
jAngDataDpt_m2 = squeeze(nanmean(jAngDataDpt_m1,6));
jTrqDatapt_m2 = squeeze(nanmean(jTrqDatapt_m1,6));
jPowDatapt_m2 = squeeze(nanmean(jPowDatapt_m1,6));
% jEgyNegDatapt_m2 = squeeze(nanmean(jEgyNegDatapt_m1,6));
% jEgyPosDatapt_m2 = squeeze(nanmean(jEgyPosDatapt_m1,6));
% comDataDpt_m2 = squeeze(nanmean(comDataDpt_m1,6));
emgDatapt_m2 = squeeze(nanmean(emgDatapt_m1,5));
% fDatapt_m2 = squeeze(nanmean(fDatapt_m1,5));

tDatabe_m2 = squeeze(nanmean(tDatabe_m1,3));
tDatape_m2 = squeeze(nanmean(tDatape_m1,4));


% Averages of within-subject standard deviations
jAngDatabt_std2b = squeeze(nanmean(jAngDatabt_std1,5));
jAngDataDbt_std2b = squeeze(nanmean(jAngDataDbt_std1,5));
jTrqDatabt_std2b = squeeze(nanmean(jTrqDatabt_std1,5));
jPowDatabt_std2b = squeeze(nanmean(jPowDatabt_std1,5));
% jEgyNegDatabt_std2b = squeeze(nanmean(jEgyNegDatabt_std1,5));
% jEgyPosDatabt_std2b = squeeze(nanmean(jEgyPosDatabt_sdt1,5));
% comDataDbt_std2b = squeeze(nanmean(comDataDbt_sdt1,5));
emgDatabt_std2b = squeeze(nanmean(emgDatabt_std1,4));

jAngDatapt_std2b = squeeze(nanmean(jAngDatapt_std1,6));
jAngDataDpt_std2b = squeeze(nanmean(jAngDataDpt_std1,6));
jTrqDatapt_std2b = squeeze(nanmean(jTrqDatapt_std1,6));
jPowDatapt_std2b = squeeze(nanmean(jPowDatapt_std1,6));
% jEgyNegDatapt_std2b = squeeze(nanmean(jEgyNegDatapt_std1,6));
% jEgyPosDatapt_std2b = squeeze(nanmean(jEgyPosDatapt_std1,6));
% comDataDpt_std2b = squeeze(nanmean(comDataDpt_std1,6));
emgDatapt_std2b = squeeze(nanmean(emgDatapt_std1,5));



% Standard deviations
jAngDatabt_std2 = squeeze(nanstd(jAngDatabt_m1,[],5));
jAngDataDbt_std2 = squeeze(nanstd(jAngDataDbt_m1,[],5));
jTrqDatabt_std2 = squeeze(nanstd(jTrqDatabt_m1,[],5));
jPowDatabt_std2 = squeeze(nanstd(jPowDatabt_m1,[],5));
% jEgyNegDatabt_std2 = squeeze(nanstd(jEgyNegDatabt_m1,[],5));
% jEgyPosDatabt_std2 = squeeze(nanstd(jEgyPosDatabt_m1,[],5));
% comDataDbt_std2 = squeeze(nanstd(comDataDbt_m1,[],5));
emgDatabt_std2 = squeeze(nanstd(emgDatabt_m1,[],4));

jAngDatapt_std2 = squeeze(nanstd(jAngDatapt_m1,[],6));
jAngDataDpt_std2 = squeeze(nanstd(jAngDataDpt_m1,[],6));
jTrqDatapt_std2 = squeeze(nanstd(jTrqDatapt_m1,[],6));
jPowDatapt_std2 = squeeze(nanstd(jPowDatapt_m1,[],6));
% jEgyNegDatapt_std2 = squeeze(nanmean(jEgyNegDatapt_m1,[],6));
% jEgyPosDatapt_std2 = squeeze(nanmean(jEgyPosDatapt_m1,[],6));
% comDataDpt_std2 = squeeze(nanstd(comDataDpt_m1,[],6));
emgDatapt_std2 = squeeze(nanstd(emgDatapt_m1,[],5));

tDatabe_std2 = squeeze(nanstd(tDatabe_m1,[],3));
tDatape_std2 = squeeze(nanstd(tDatape_m1,[],4));

%% Project average resmapled time-series back on time axis
% Using repetition average durations of the sequences
% This is non-cumulative, so the time axis resets at every phase
% (i.e. starts at 0)

% Convert event times to matching sequence durations
tDatabt_m2 = zeros(7,4);
tDatabt_m2(1,:) = 0.15; % TOR-PertEnd
tDatabt_m2(2,:) = tDatabe_m2(3,:) - 0.15; % PertEnd-HSR
tDatabt_m2(3:6,:) = tDatabe_m2(3:6,:); % Remaining sequence durations
tDatabt_m2(7,:) = sum(tDatabe_m2(3:6,:),1); % TOR-TOR2
tDatabt_m2 = repmat( [(0:nsmpl-1)/nsmpl].' , [1 4 7]) .* permute( repmat( tDatabt_m2 , [1 1 nsmpl]) , [3 2 1]);

tDatapt_m2 = zeros(7,8,4);
tDatapt_m2(1,:,:) = 0.15; % TOR-PertEnd
tDatapt_m2(2,:,:) = tDatape_m2(3,:,:) - 0.15; % PertEnd-HSR
tDatapt_m2(3:6,:,:) = tDatape_m2(3:6,:,:); % Remaining sequence durations
tDatapt_m2(7,:,:) = sum(tDatape_m2(3:6,:,:),1); % TOR-TOR2
tDatapt_m2 = repmat( [(0:nsmpl-1)/nsmpl].' , [1 8 4 7]) .* permute( repmat( tDatapt_m2, [1 1 1 nsmpl]) , [4 2 3 1]);

