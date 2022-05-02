%% Plotting


%%
% #######################
% ## Grouped subjects  ##
% #######################


%% Time

close all; clc;

chanpTtl = {'SlowML','FastML','SlowAP','FastAP'};
altord = [4:-1:1 6:9];

for iset = 1:4
    figure; 
    
    h1 = subplot(121); hold on;
    
    % BASE
    errorbar(5,tDatabe_m2(3,iset),tDatabe_std2(3,iset),'color',colb);
    plot(5,tDatabe_m2(3,iset),'o','color',colb,'markerfacecolor',colb);
    
    errorbar(5,tDatabe_m2(4,iset),tDatabe_std2(4,iset),'color',colb);
    plot(5,tDatabe_m2(4,iset),'o','color',colb);
    
    % PERT
    for ipert = 1:8
        
        errorbar(altord(ipert),tDatape_m2(3,ipert,iset),tDatape_std2(3,ipert,iset),'color',col(ipert,:));
        plot(altord(ipert),tDatape_m2(3,ipert,iset),'o','color',col(ipert,:),'markerfacecolor',col(ipert,:));
        
        errorbar(altord(ipert),tDatape_m2(4,ipert,iset),tDatape_std2(4,ipert,iset),'color',col(ipert,:));
        plot(altord(ipert),tDatape_m2(4,ipert,iset),'o','color',col(ipert,:));
        
    end
    axis tight;
    axis square;
    
    h2 = subplot(122); hold on;
    
    % BASE
    errorbar(5,tDatabe_m2(5,iset),tDatabe_std2(5,iset),'color',colb);
    plot(5,tDatabe_m2(5,iset),'o','color',colb,'markerfacecolor',colb);
    
    errorbar(5,tDatabe_m2(6,iset),tDatabe_std2(6,iset),'color',colb);
    plot(5,tDatabe_m2(6,iset),'o','color',colb);
    
    % PERT
    for ipert = 1:8
        
        errorbar(altord(ipert),tDatape_m2(5,ipert,iset),tDatape_std2(5,ipert,iset),'color',col(ipert,:));
        plot(altord(ipert),tDatape_m2(5,ipert,iset),'o','color',col(ipert,:),'markerfacecolor',col(ipert,:));
        
        errorbar(altord(ipert),tDatape_m2(6,ipert,iset),tDatape_std2(6,ipert,iset),'color',col(ipert,:));
        plot(altord(ipert),tDatape_m2(6,ipert,iset),'o','color',col(ipert,:));
        
    end
    axis tight;
    axis square;
    
    % Ylim
    ymin = min(cellfun(@min,get([h1 h2],'ylim')));
    ymax = max(cellfun(@max,get([h1 h2],'ylim')));
    set([h1 h2],'xlim',[0.5 9.5],'ylim',[ymin ymax]);
    
    % Ticks & labels
    set([h1 h2],'xtick',[1:2:9],'xticklabel',[-0.16 -0.08 0 0.08 0.16]);
    set(h2,'ytick',[])
    xlabel(h1,'left/backward <  > right/forward');
    xlabel(h2,'left/backward <  > right/forward');
    ylabel(h1,'Duration [s]');
    
    % Titles
    title(h1,'TOR-HSR & HSR-TOL');
    title(h2,'TOL-HSL & HSL-TOR');
    suptitle(['Pert:' chanpTtl{iset}]);
    
    % Print figure
    fullfilename = [pwd '\Figs\Times\' chanpTtl{iset} '_COMD'];
    print(gcf,fullfilename,'-r300','-dpng')
    
end


%% COM velocity

close all;
clear hplt;
chanpTtl = {'SlowML','FastML','SlowAP','FastAP'};

seqs = [1 2 4:6];

for iset = 1:4
    figure;
    for idim = 1:3
        for iseq = seqs

            switch iseq
                case 1
                    hplt(1+5*(idim-1)) = subplot(3,6,1+6*(idim-1)); hold on;
                case 2
                    hplt(2+5*(idim-1)) = subplot(3,6,2+6*(idim-1)); hold on;
                case 4
                    hplt(3+5*(idim-1)) = subplot(3,6,3+6*(idim-1)); hold on;
                case 5
                    hplt(4+5*(idim-1)) = subplot(3,6,[4 5]+6*(idim-1)); hold on;
                case 6
                    hplt(5+5*(idim-1)) = subplot(3,6,6+6*(idim-1)); hold on;
            end
                
            % BASE
            plot(comDataDbt_m2(:,end,idim,iset,iseq),'-','color','k','linewidth',2);

            % PERT
            for ipert = 1:8
                plot(comDataDpt_m2(:,end,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1.5);
            end

            axis tight;
            
        end
        
        % xy limits
        ymin = min(cellfun(@min,get(hplt((1:5)+5*(idim-1)),'ylim')));
        ymax = max(cellfun(@max,get(hplt((1:5)+5*(idim-1)),'ylim')));
        set(hplt((1:5)+5*(idim-1)),'xlim',[1 50],'ylim',[ymin ymax]);

    end

    % xy labels
    set(hplt,'xtick',[]);
    xlabel(hplt(11),'TOR(PStart)-PEnd');
    xlabel(hplt(12),'PEnd-HSR');
    xlabel(hplt(13),'HSR-TOL');
    xlabel(hplt(14),'TOL-HSL');
    xlabel(hplt(15),'HSL-TOR');
    ylabel(hplt(1),'ML (m/s)');
    ylabel(hplt(6),'AP  (m/s)');
    ylabel(hplt(11),'VT  (m/s)');
    set(hplt([2:5 7:10 12:15]),'ytick',[],'ycolor',[1 1 1]);

    % Title
    suptitle(['Pert:' chanpTtl{iset}]);

    % Print figure
    fullfilename = [pwd '\Figs\COMD\' chanpTtl{iset} '_COMD'];
    print(gcf,fullfilename,'-r300','-dpng')

end

%% EMG signals (as resampled series)
% 
%  'TAL'    'GML'    'RFL'    'BFL'    'GAML'    'ALL'    'TAR'    'GMR'    'RFR'    'BFR'    'GAMR'    'ALR'

close all;
clear hplt;

chanjTtl = {'TibAnt','GasMed','RecFem','BicFem','GluMed','AddLon'};
chanpTtl = {'SlowML','FastML','SlowAP','FastAP'};
seqs = [1 2 4:6];

for iset = 1:4

    for ichan = 1:6
        figure;

        for iside = 0:1
            for iseq = seqs

                switch iseq
                    case 1
                        hplt(1+5*iside) = subplot(2,6,1+6*iside); hold on;
                    case 2
                        hplt(2+5*iside) = subplot(2,6,2+6*iside); hold on;
                    case 4
                        hplt(3+5*iside) = subplot(2,6,3+6*iside); hold on;
                    case 5
                        hplt(4+5*iside) = subplot(2,6,[4 5]+6*iside); hold on;
                    case 6
                        hplt(5+5*iside) = subplot(2,6,6+6*iside); hold on;
                end

                % BASE
                plot( emgDatabt_m2(:,ichan+6*iside,iset,iseq),'color','k','linewidth',2);
                fill( [1:size(emgDatabt_m2,1) size(emgDatabt_m2,1):-1:1],...
                    [emgDatabt_m2(:,ichan+6*iside,iset,iseq)+emgDatabt_std2(:,ichan+6*iside,iset,iseq) ; ...
                    emgDatabt_m2(end:-1:1,ichan+6*iside,iset,iseq)-emgDatabt_std2(end:-1:1,ichan+6*iside,iset,iseq)],...
                    colb,'edgecolor','none','facealpha',0.4);
                
                % PERT
                for ipert = 1:8
                    plot( emgDatapt_m2(:,ichan+6*iside,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1.5);
                end

                axis tight;

            end
        end
        
        % xy limits
        ymin = min(cellfun(@min,get(hplt,'ylim')));
        ymax = max(cellfun(@max,get(hplt,'ylim')));
        set(hplt,'xlim',[1 50],'ylim',[ymin ymax]);

        % xy labels
        set(hplt,'xtick',[]);
        xlabel(hplt(6),'TOR(PStart)-PEnd');
        xlabel(hplt(7),'PEnd-HSR');
        xlabel(hplt(8),'HSR-TOL');
        xlabel(hplt(9),'TOL-HSL');
        xlabel(hplt(10),'HSL-TOR');
        ylabel(hplt(1),'LEFT');
        ylabel(hplt(6),'RIGHT');
        set(hplt([2:5 7:10]),'ytick',[],'ycolor',[1 1 1]);

        % Title
        suptitle(['Pert:' chanpTtl{iset} ' | Muscle:' chanjTtl{ichan}]);

        % Print figure
        fullfilename = [pwd '\Figs\EMG\' chanpTtl{iset} '_' chanjTtl{ichan}];
        print(gcf,fullfilename,'-r300','-dpng')
        
    end

end

%% Joint angles
%
% HipAP pos is (ante)flexion, neg is extension
% HipML pos is adduction, neg is abduction
% KneAP pos is extension, neg is flexion
% AnkAP pos is dorsiflexion, neg is plantarflexion
% AnkML pos is eversion, neg is inversion
% 
% Note that changes in hip joint angle can also occur due to pelvis rotation!

close all;
clear hplt;

% iset = 4;
sets = 1;
for iset = sets

    
seqs = [1 2 4:6];
chandimL = [4 4 5 6 6 ; 1 2 1 1 2];
chandimR = [7 7 8 9 9 ; 1 2 1 1 2];
chanjTtl = {'HipAP','HipML','KneAP','AnkAP','AnkML'};
chanpTtl = {'SlowML','FastML','SlowAP','FastAP'};
chanXlbl = {{'Ext <> Flex','Rad'},{'Abd <> Add','Rad'},{'Flex <> Ext','Rad'},{'PFlex <> DFlex','Rad'},{'Inver <> Ever','Rad'}};

for iplt = 1:size(chandimL,2)
    
    figure('name',[chanpTtl{iset} '|' chanjTtl{iplt}],...
        'color',[1 1 1],...
        'position',[400 150 1080 720]);
    
    for iside = 0:1
        
        if iside == 0
            ichan = chandimL(1,iplt);
            idim = chandimL(2,iplt);
        else
            ichan = chandimR(1,iplt);
            idim = chandimR(2,iplt);
        end
        
        for iseq = seqs

            switch iseq
                case 1
                    hplt(1+5*iside) = subplot(2,6,1+6*iside); hold on;
                case 2
                    hplt(2+5*iside) = subplot(2,6,2+6*iside); hold on;
                case 4
                    hplt(3+5*iside) = subplot(2,6,3+6*iside); hold on;
                case 5
                    hplt(4+5*iside) = subplot(2,6,[4 5]+6*iside); hold on;
                case 6
                    hplt(5+5*iside) = subplot(2,6,6+6*iside); hold on;
            end

            % BASE
            plot(jAngDatabt_m2(:,ichan,idim,iset,iseq),'-','linewidth',2,'color','k');
%             fill([1:size(jAngDatabt_m2,1) size(jAngDatabt_m2,1):-1:1],...
%                 [jAngDatabt_m2(:,ichan,idim,iset,iseq)+jAngDatabt_std2(:,ichan,idim,iset,iseq) ;...
%                 jAngDatabt_m2(end:-1:1,ichan,idim,iset,iseq)-jAngDatabt_std2(end:-1:1,ichan,idim,iset,iseq)],...
%                 colb,'facealpha',0.4,'edgecolor','none');
                
%                 plot(jAngDataDbt_m2(:,ichan,idim,iset,iseq),'-','linewidth',2,'color','k'); % VELOCITIES
%                 plot(jAngDatabt_m2(:,2,3,1,iseq),'-','linewidth',2,'color','k')
                
            % PERT
            for ipert = 1:8
                plot(jAngDatapt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1.5);
    %             fill([1:size(jAngDatapt_m2) size(jAngDatapt_m2):-1:1],...
    %                 [jAngDatapt_m2(:,ichan,idim,ipert,iset,iseq)+jAngDatapt_std2(:,ichan,idim,ipert,iset,iseq) ;...
    %                 jAngDatapt_m2(end:-1:1,ichan,idim,ipert,iset,iseq)-jAngDatapt_std2(end:-1:1,ichan,idim,ipert,iset,iseq)],...
    %                 col(ipert,:),'facealpha',0.1,'edgecolor','none');
    
%                     plot(jAngDataDpt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1.5); % VELOCITIES
%                     plot(jAngDatapt_m2(:,2,3,ipert,1,iseq),'-','color',col(ipert,:),'linewidth',1.5);
            end
            
            axis tight;
            
        end
    end

    % xy limits
    ymin = min(cellfun(@min,get(hplt,'ylim')));
    ymax = max(cellfun(@max,get(hplt,'ylim')));
    set(hplt,'xlim',[1 50],'ylim',[ymin ymax]);
    
    % xy labels
    set(hplt,'xtick',[]);
    xlabel(hplt(6),'TOR(PStart)-PEnd');
    xlabel(hplt(7),'PEnd-HSR');
    xlabel(hplt(8),'HSR-TOL');
    xlabel(hplt(9),'TOL-HSL');
    xlabel(hplt(10),'HSL-TOR');
    ylabel(hplt(1),chanXlbl{iplt});
    ylabel(hplt(6),chanXlbl{iplt});
    set(hplt([2:5 7:10]),'ytick',[],'ycolor',[1 1 1]);
    
    % Title
    suptitle(['Pert:' chanpTtl{iset} ' | Joint:' chanjTtl{iplt}]);

    % Print figure
%     fullfilename = [pwd '\Figs\JAng\' chanpTtl{iset} '_' chanjTtl{iplt}];
%     print(gcf,fullfilename,'-r300','-dpng')

end

end


%% Joint angles, on averaged time axis
%
% HipAP pos is (ante)flexion, neg is extension
% HipML pos is adduction, neg is abduction
% KneAP pos is extension, neg is flexion
% AnkAP pos is dorsiflexion, neg is plantarflexion
% AnkML pos is eversion, neg is inversion
% 
% Note that changes in hip joint angle can also occur due to pelvis rotation!

close all;
clear hplt;

sets = 4;
% sets = 1;
for iset = sets

    
seqs = [1 2 4:6];
chandimL = [4 4 5 6 6 ; 1 2 1 1 2];
chandimR = [7 7 8 9 9 ; 1 2 1 1 2];
chanjTtl = {'HipAP','HipML','KneAP','AnkAP','AnkML'};
chanpTtl = {'SlowML','FastML','SlowAP','FastAP'};
chanXlbl = {{'Ext <> Flex','Rad'},{'Abd <> Add','Rad'},{'Flex <> Ext','Rad'},{'PFlex <> DFlex','Rad'},{'Inver <> Ever','Rad'}};

for iplt = 1:size(chandimL,2)
    
    figure('name',[chanpTtl{iset} '|' chanjTtl{iplt}],...
        'color',[1 1 1],...
        'position',[400 150 1080 720]);
    
    for iside = 0:1
        
        if iside == 0
            ichan = chandimL(1,iplt);
            idim = chandimL(2,iplt);
        else
            ichan = chandimR(1,iplt);
            idim = chandimR(2,iplt);
        end
        
        for iseq = seqs

            switch iseq
                case 1
                    hplt(1+5*iside) = subplot(2,6,1+6*iside); hold on;
                case 2
                    hplt(2+5*iside) = subplot(2,6,2+6*iside); hold on;
                case 4
                    hplt(3+5*iside) = subplot(2,6,3+6*iside); hold on;
                case 5
                    hplt(4+5*iside) = subplot(2,6,[4 5]+6*iside); hold on;
                case 6
                    hplt(5+5*iside) = subplot(2,6,6+6*iside); hold on;
            end

            % BASE
            plot(tDatabt_m2(:,iset,iseq),jAngDatabt_m2(:,ichan,idim,iset,iseq),'-','linewidth',2,'color','k');
%             fill([1:size(jAngDatabt_m2,1) size(jAngDatabt_m2,1):-1:1],...
%                 [jAngDatabt_m2(:,ichan,idim,iset,iseq)+jAngDatabt_std2(:,ichan,idim,iset,iseq) ;...
%                 jAngDatabt_m2(end:-1:1,ichan,idim,iset,iseq)-jAngDatabt_std2(end:-1:1,ichan,idim,iset,iseq)],...
%                 colb,'facealpha',0.4,'edgecolor','none');
                
%                 plot(jAngDataDbt_m2(:,ichan,idim,iset,iseq),'-','linewidth',2,'color','k'); % VELOCITIES
%                 plot(jAngDatabt_m2(:,2,3,1,iseq),'-','linewidth',2,'color','k')

            % PERT
            for ipert = 1:8
                plot(tDatapt_m2(:,ipert,iset,iseq),jAngDatapt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1.5);
    %             fill([1:size(jAngDatapt_m2) size(jAngDatapt_m2):-1:1],...
    %                 [jAngDatapt_m2(:,ichan,idim,ipert,iset,iseq)+jAngDatapt_std2(:,ichan,idim,ipert,iset,iseq) ;...
    %                 jAngDatapt_m2(end:-1:1,ichan,idim,ipert,iset,iseq)-jAngDatapt_std2(end:-1:1,ichan,idim,ipert,iset,iseq)],...
    %                 col(ipert,:),'facealpha',0.1,'edgecolor','none');
    
%                     plot(jAngDataDpt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1.5); % VELOCITIES
%                     plot(jAngDatapt_m2(:,2,3,ipert,1,iseq),'-','color',col(ipert,:),'linewidth',1.5);
            end
            
            axis tight;
            
        end
    end

    % xy limits
    ymin = min(cellfun(@min,get(hplt,'ylim')));
    ymax = max(cellfun(@max,get(hplt,'ylim')));
%     set(hplt,'xlim',[1 50],'ylim',[ymin ymax]);
    set(hplt,'ylim',[ymin ymax]);
    
    % xy labels
%     set(hplt,'xtick',[]);
    xlabel(hplt(6),'TOR(PStart)-PEnd');
    xlabel(hplt(7),'PEnd-HSR');
    xlabel(hplt(8),'HSR-TOL');
    xlabel(hplt(9),'TOL-HSL');
    xlabel(hplt(10),'HSL-TOR');
    ylabel(hplt(1),chanXlbl{iplt});
    ylabel(hplt(6),chanXlbl{iplt});
    set(hplt([2:5 7:10]),'ytick',[],'ycolor',[1 1 1]);
    
    % Title
    suptitle(['Pert:' chanpTtl{iset} ' | Joint:' chanjTtl{iplt}]);

    % Print figure
    fullfilename = [pwd '\Figs\JAngTimeAx\' chanpTtl{iset} '_' chanjTtl{iplt}];
    print(gcf,fullfilename,'-r300','-dpng')

end

end

%% Joint angles, torques, and power on averaged time axis, double column
%
% HipAP pos is (ante)flexion, neg is extension
% HipML pos is adduction, neg is abduction
% KneAP pos is extension, neg is flexion
% AnkAP pos is dorsiflexion, neg is plantarflexion
% AnkML pos is eversion, neg is inversion
% 
% Note that changes in hip joint angle can also occur due to pelvis rotation!

close all;
clear hplt;

% sets = 4;
sets = 1:4;
for iset = sets

seqs = [1 2 4:6];
chandimL = [4 4 5 6 6 ; 1 2 1 1 2];
chandimR = [7 7 8 9 9 ; 1 2 1 1 2];
% chandimL = [4 4 5 6 6 ; 3 3 3 3 3];
% chandimR = [7 7 8 9 9 ; 3 3 3 3 3];
chanjTtl = {'HipAP','HipML','KneAP','AnkAP','AnkML'};
chanpTtl = {'SlowML','FastML','SlowAP','FastAP'};
% chanXlblA = {{'Ext <> Flex','Rad'},{'Abd <> Add','Rad'},{'Flex <> Ext','Rad'},{'PFlex <> DFlex','Rad'},{'Inver <> Ever','Rad'}};
% chanXlblT = {{'Ext <> Flex','Nm/mgl'},{'Abd <> Add','Nm/mgl'},{'Flex <> Ext','Nm/mgl'},{'PFlex <> DFlex','Nm/mgl'},{'Inver <> Ever','Nm/mgl'}};
% chanYlblA = {{'Angle','Ext <> Flex'},{'Angle','Abd <> Add'},{'Angle','Flex <> Ext'},{'Angle','PFlex <> DFlex'},{'Angle','Inver <> Ever'}};

jTtl = {'Hip flexion-extension','Hip abduction-adduction','Knee flexion-extension','Ankle plantarflexion-dorsiflexion','Ankle inversion-eversion'};
pTtl = {'Slow walking, ML perturbations','Normal walking, ML perturbations','Slow walking, AP perturbations','Normal walking, AP perturbations'};
chanYlblA = {'Ext  |  Angle  |  Flex','Abd  |  Angle  |  Add','Flex  |  Angle  |  Ext','PFlex  |  Angle  |  DFlex','Inver  |  Angle  |  Ever'};
chanYlblT = {'Ext  |  Torque  |  Flex','Abd  |  Torque  |  Add','Flex  |  Torque  |  Ext','PFlex  |  Torque  |  DFlex','Inver  |  Torque  |  Ever'};


for iplt = 1:size(chandimL,2)
    
%     figure('name',[chanpTtl{iset} '|' chanjTtl{iplt}],...
%         'color',[1 1 1],...
%         'position',[220 126 1360 850]);

%     figure('name',[chanpTtl{iset} '|' chanjTtl{iplt}],...
%         'color',[1 1 1],...
%         'position',[533 404 850 529]);  % Without suptitle
    
    figure('name',[chanpTtl{iset} '|' chanjTtl{iplt}],...
        'color',[1 1 1],...
        'position',[533 404 850 529+0.1*529]);  % With suptitle
    
    for itype = 0:2 % angle (0), torque (1), and power (2)
    
        for iside = 0:1
            if iside == 0
                ichan = chandimL(1,iplt);
                idim = chandimL(2,iplt);
            else
                ichan = chandimR(1,iplt);
                idim = chandimR(2,iplt);
            end

            for iseq = seqs

                switch iseq
                    case 1
                        hplt(1+5*iside+10*itype) = subplot(3,12,1+6*iside+12*itype); hold on;
                    case 2
                        hplt(2+5*iside+10*itype) = subplot(3,12,2+6*iside+12*itype); hold on;
                    case 4
                        hplt(3+5*iside+10*itype) = subplot(3,12,3+6*iside+12*itype); hold on;
                    case 5
                        hplt(4+5*iside+10*itype) = subplot(3,12,[4 5]+6*iside+12*itype); hold on;
                    case 6
                        hplt(5+5*iside+10*itype) = subplot(3,12,6+6*iside+12*itype); hold on;
                end

                if itype == 0 % JOINT ANGLES
                    
                    % BASE
                    fill([tDatabt_m2(:,iset,iseq) ; tDatabt_m2(end:-1:1,iset,iseq)],...
                        [jAngDatabt_m2(:,ichan,idim,iset,iseq)+jAngDatabt_std2b(:,ichan,idim,iset,iseq) ;...
                        jAngDatabt_m2(end:-1:1,ichan,idim,iset,iseq)-jAngDatabt_std2b(end:-1:1,ichan,idim,iset,iseq)],...
                        colb,'edgecolor','none');
                    plot(tDatabt_m2(:,iset,iseq),jAngDatabt_m2(:,ichan,idim,iset,iseq),'-','linewidth',1,'color','k');

        %                 plot(jAngDataDbt_m2(:,ichan,idim,iset,iseq),'-','linewidth',2,'color','k'); % VELOCITIES
        %                 plot(jAngDatabt_m2(:,2,3,1,iseq),'-','linewidth',2,'color','k')

                    % PERT
                    for ipert = 1:8
%                         if any(ipert == [4 8])
%                             fill([tDatapt_m2(:,ipert,iset,iseq) ; tDatapt_m2(end:-1:1,ipert,iset,iseq)],...
%                                 [jAngDatapt_m2(:,ichan,idim,ipert,iset,iseq)+jAngDatapt_std2b(:,ichan,idim,ipert,iset,iseq) ;...
%                                 jAngDatapt_m2(end:-1:1,ichan,idim,ipert,iset,iseq)-jAngDatapt_std2b(end:-1:1,ichan,idim,ipert,iset,iseq)],...
%                                 col(ipert,:),'facealpha',0.1,'edgecolor','none');
%                         end

                        plot(tDatapt_m2(:,ipert,iset,iseq),jAngDatapt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1);
    %                     plot(jAngDataDpt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1.5); % VELOCITIES
    %                     plot(jAngDatapt_m2(:,2,3,ipert,1,iseq),'-','color',col(ipert,:),'linewidth',1.5);
                    end

                
                elseif itype == 1 % JOINT TORQUES
                    
                    % BASE
                    fill([tDatabt_m2(:,iset,iseq) ; tDatabt_m2(end:-1:1,iset,iseq)],...
                        [jTrqDatabt_m2(:,ichan,idim,iset,iseq)+jTrqDatabt_std2b(:,ichan,idim,iset,iseq) ;...
                        jTrqDatabt_m2(end:-1:1,ichan,idim,iset,iseq)-jTrqDatabt_std2b(end:-1:1,ichan,idim,iset,iseq)],...
                        colb,'edgecolor','none');
                    plot(tDatabt_m2(:,iset,iseq),jTrqDatabt_m2(:,ichan,idim,iset,iseq),'-','linewidth',1,'color','k');

                    % PERT
                    if any(iset == [3 4])
                        setperts = 1:8;
                    else
                        setperts = 5:8; % No inward
                    end
                    for ipert = setperts
%                         fill([tDatapt_m2(:,ipert,iset,iseq) ; tDatapt_m2(end:-1:1,ipert,iset,iseq)],...
%                             [jTrqDatapt_m2(:,ichan,idim,ipert,iset,iseq)+jTrqDatapt_std2(:,ichan,idim,ipert,iset,iseq) ;...
%                             jTrqDatapt_m2(end:-1:1,ichan,idim,ipert,iset,iseq)-jTrqDatapt_std2(end:-1:1,ichan,idim,ipert,iset,iseq)],...
%                             col(ipert,:),'facealpha',0.1,'edgecolor','none');
                        plot(tDatapt_m2(:,ipert,iset,iseq),jTrqDatapt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1);
                    end
                    
                else % JOINT POWER
                    
                    % BASE
                    fill([tDatabt_m2(:,iset,iseq) ; tDatabt_m2(end:-1:1,iset,iseq)],...
                        [jPowDatabt_m2(:,ichan,idim,iset,iseq)+jPowDatabt_std2b(:,ichan,idim,iset,iseq) ;...
                        jPowDatabt_m2(end:-1:1,ichan,idim,iset,iseq)-jPowDatabt_std2b(end:-1:1,ichan,idim,iset,iseq)],...
                        colb,'edgecolor','none');
                    plot(tDatabt_m2(:,iset,iseq),jPowDatabt_m2(:,ichan,idim,iset,iseq),'-','linewidth',1,'color','k');

                    % PERT
                    if any(iset == [3 4])
                        setperts = 1:8;
                    else
                        setperts = 5:8; % No inward
                    end
                    for ipert = setperts
%                         fill([tDatapt_m2(:,ipert,iset,iseq) ; tDatapt_m2(end:-1:1,ipert,iset,iseq)],...
%                             [jTrqDatapt_m2(:,ichan,idim,ipert,iset,iseq)+jTrqDatapt_std2(:,ichan,idim,ipert,iset,iseq) ;...
%                             jTrqDatapt_m2(end:-1:1,ichan,idim,ipert,iset,iseq)-jTrqDatapt_std2(end:-1:1,ichan,idim,ipert,iset,iseq)],...
%                             col(ipert,:),'facealpha',0.1,'edgecolor','none');
                        plot(tDatapt_m2(:,ipert,iset,iseq),jPowDatapt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1);
                    end
                    
                end

                % Tight axes
                axis tight;
                set(gca,'color','none')
                
                % Plot dashed zero line (after tight)
                maxx = max([tDatabt_m2(:,iset,iseq) ; reshape(tDatapt_m2(:,:,iset,iseq),[400,1])]);
                plot([0 maxx],[0 0],'--k');

            end
        end
    
    end

    % xy limits
    ymin = min(cellfun(@min,get(hplt(1:10),'ylim')));
    ymax = max(cellfun(@max,get(hplt(1:10),'ylim')));
    set(hplt(1:10),'ylim',[ymin ymax]);
    
    ymin = min(cellfun(@min,get(hplt(11:20),'ylim')));
    ymax = max(cellfun(@max,get(hplt(11:20),'ylim')));
    set(hplt(11:20),'ylim',[ymin ymax]);
    
    ymin = min(cellfun(@min,get(hplt(21:30),'ylim')));
    ymax = max(cellfun(@max,get(hplt(21:30),'ylim')));
    set(hplt(21:30),'ylim',[ymin ymax]);
    
    % Suptitle
    httl = suptitle(sprintf([pTtl{iset} '\n' jTtl{iplt}]));
    set(httl,'fontname','arial','fontsize',12);
    
    % xy labels
    set(hplt,'fontname','Arial');
    
    set(hplt,'xtick',[]);
    ytix = get(hplt(1),'ytick');
    set(hplt(1),'fontsize',8,'ytick',ytix);
%     set(hplt(1),'yticklabel',cellstr(num2str(ytix(:)))); % Convert to string to prevent double text objects when exporting
    set(hplt([2:10]),'yticklabel',[]);
    set(hplt([2:5 7:10]),'ycolor',[1 1 1]);
    
    set(hplt,'xtick',[]);
    ytix = get(hplt(11),'ytick');
    set(hplt(11),'fontsize',8,'ytick',ytix);
    set(hplt([12:20]),'yticklabel',[]);
    set(hplt([12:15 17:20]),'ycolor',[1 1 1]);
    
    set(hplt,'xtick',[]);
    ytix = get(hplt(21),'ytick');
    set(hplt(21),'fontsize',8,'ytick',ytix);
    set(hplt([22:30]),'yticklabel',[]);
    set(hplt([22:25 27:30]),'ycolor',[1 1 1]);
    
    ylabel(hplt(1),chanYlblA{iplt});
        ylab = get(hplt(1),'ylabel');
        ylabp = get(ylab,'position');
        set(ylab,'position',[-0.12 ylabp(2:3)],'fontsize',8);
    ylabel(hplt(11),chanYlblT{iplt});
        ylab = get(hplt(11),'ylabel');
        ylabp = get(ylab,'position');
        set(ylab,'position',[-0.12 ylabp(2:3)],'fontsize',8);
    ylabel(hplt(21),'Power');
        ylab = get(hplt(21),'ylabel');
        ylabp = get(ylab,'position');
        set(ylab,'position',[-0.12 ylabp(2:3)],'fontsize',8);
        
    % Print figures    
    fullfilename = [pwd '\Figs\ExportFolder\' chanpTtl{iset} '_' chanjTtl{iplt}];
%         
%     % Print figure as png
%     print(gcf,fullfilename,'-r300','-dpng')
% 
    % Print figure as eps
    printeps(get(gcf,'number'),fullfilename);

end
% return;
end

%% EMG, double column

%  'TAL'    'GML'    'RFL'    'BFL'    'GAML'    'ALL'    'TAR'    'GMR'    'RFR'    'BFR'    'GAMR'    'ALR'

close all;
clear hplt;

chanmTtl = {'TA','GM','RF','BF','GLM','AL'};
chanpTtl = {'SlowML','FastML','SlowAP','FastAP'};
pTtl = {'Slow walking, ML perturbations','Normal walking, ML perturbations','Slow walking, AP perturbations','Normal walking, AP perturbations'};
seqs = [1 2 4:6];

for iset = 1:4

    clear hplt
    
    figure('name',chanpTtl{iset},...
    'color',[1 1 1],...
    'position',[533 404 850 529+0.1*529]);  % With suptitle
    
    for ichan = 1:6

        for iside = 0:1
            for iseq = seqs

                switch iseq
                    case 1
                        hplt(1+5*iside+10*(ichan-1)) = subplot(6,12,1+6*iside+12*(ichan-1)); hold on;
                    case 2
                        hplt(2+5*iside+10*(ichan-1)) = subplot(6,12,2+6*iside+12*(ichan-1)); hold on;
                    case 4
                        hplt(3+5*iside+10*(ichan-1)) = subplot(6,12,3+6*iside+12*(ichan-1)); hold on;
                    case 5
                        hplt(4+5*iside+10*(ichan-1)) = subplot(6,12,[4 5]+6*iside+12*(ichan-1)); hold on;
                    case 6
                        hplt(5+5*iside+10*(ichan-1)) = subplot(6,12,6+6*iside+12*(ichan-1)); hold on;
                end

                % BASE
                fill( [tDatabt_m2(:,iset,iseq) ; tDatabt_m2(end:-1:1,iset,iseq)],...
                    [emgDatabt_m2(:,ichan+6*iside,iset,iseq)+emgDatabt_std2(:,ichan+6*iside,iset,iseq) ; ...
                    emgDatabt_m2(end:-1:1,ichan+6*iside,iset,iseq)-emgDatabt_std2(end:-1:1,ichan+6*iside,iset,iseq)],...
                    colb,'edgecolor','none');
                plot( tDatabt_m2(:,iset,iseq) , emgDatabt_m2(:,ichan+6*iside,iset,iseq),'color','k','linewidth',1);
                
                % PERT
                for ipert = 1:8
                    plot( tDatapt_m2(:,ipert,iset,iseq) , emgDatapt_m2(:,ichan+6*iside,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1);
                end

                % Tight axes
                axis tight;
                set(gca,'color','none')

            end
        end
    end

        % xy limits
        for iplt = 1:6
            
            ymin = min(cellfun(@min,get(hplt((1:10)+10*(iplt-1)),'ylim'))); % Row 1
            ymax = max(cellfun(@max,get(hplt((1:10)+10*(iplt-1)),'ylim')));
            set(hplt((1:10)+10*(iplt-1)),'ylim',[ymin ymax]);

%             ymin = min(cellfun(@min,get(hplt((6:10)+10*(iplt-1)),'ylim')));
%             ymax = max(cellfun(@max,get(hplt((6:10)+10*(iplt-1)),'ylim')));
%             set(hplt((6:10)+10*(iplt-1)),'ylim',[ymin ymax]);
            
        end
        
        % Suptitle
        httl = suptitle(pTtl(iset));
        set(httl,'fontname','arial','fontsize',12);

        % xy labels
        set(hplt,'fontname','Arial');

        for iplt = 1:6
            
            % xy labels
            set(hplt,'xtick',[]);
            ytix = get(hplt(1+10*(iplt-1)),'ytick');
            set(hplt(1+10*(iplt-1)),'fontsize',8,'ytick',ytix);
%             ytix = get(hplt(6+10*(iplt-1)),'ytick');
%             set(hplt(6+10*(iplt-1)),'fontsize',8,'ytick',ytix);
            set(hplt([2:5 6:10]+10*(iplt-1)),'yticklabel',[]);
            set(hplt([2:5 7:10]+10*(iplt-1)),'ycolor',[1 1 1]);
        
            % Ylabel
            ylabel(hplt(1+10*(iplt-1)),chanmTtl{iplt});
                ylab = get(hplt(1+10*(iplt-1)),'ylabel');
                ylabp = get(ylab,'position');
                set(ylab,'position',[-0.12 ylabp(2:3)],'fontsize',8);
        end

        % Print figures    
        fullfilename = [pwd '\Figs\ExportFolder\EMG_' chanpTtl{iset}];
            
        % Print figure as png
        print(gcf,fullfilename,'-r300','-dpng')
    
        % Print figure as eps
        printeps(get(gcf,'number'),fullfilename);
        
% return;
end


%% Joint torques
%
% HipAP pos is (ante)flexion, neg is extension
% HipML pos is adduction, neg is abduction
% KneAP pos is extension, neg is flexion
% AnkAP pos is dorsiflexion, neg is plantarflexion
% AnkML pos is eversion, neg is inversion
% 
% Note that changes in hip joint angle can also occur due to pelvis rotation!

close all;
clear hplt;

sets = 4;%1:4;
for iset = sets

    
seqs = [1 2 4:6];
chandimL = [4 4 5 6 6 ; 1 2 1 1 2];
chandimR = [7 7 8 9 9 ; 1 2 1 1 2];
chanjTtl = {'HipAP','HipML','KneAP','AnkAP','AnkML'};
chanpTtl = {'SlowML','FastML','SlowAP','FastAP'};
chanXlbl = {{'Ext <> Flex','Nm/mgl'},{'Abd <> Add','Nm/mgl'},{'Flex <> Ext','Nm/mgl'},{'PFlex <> DFlex','Nm/mgl'},{'Inver <> Ever','Nm/mgl'}};

for iplt = 1:size(chandimL,2)
    
    figure('name',[chanpTtl{iset} '|' chanjTtl{iplt}],...
        'color',[1 1 1],...
        'position',[400 150 1080 720]);
    
    for iside = 0:1
        
        if iside == 0
            ichan = chandimL(1,iplt);
            idim = chandimL(2,iplt);
        else
            ichan = chandimR(1,iplt);
            idim = chandimR(2,iplt);
        end
        
        for iseq = seqs

            switch iseq
                case 1
                    hplt(1+5*iside) = subplot(2,6,1+6*iside); hold on;
                case 2
                    hplt(2+5*iside) = subplot(2,6,2+6*iside); hold on;
                case 4
                    hplt(3+5*iside) = subplot(2,6,3+6*iside); hold on;
                case 5
                    hplt(4+5*iside) = subplot(2,6,[4 5]+6*iside); hold on;
                case 6
                    hplt(5+5*iside) = subplot(2,6,6+6*iside); hold on;
            end

            % BASE
            plot(jTrqDatabt_m2(:,ichan,idim,iset,iseq),'-','linewidth',2,'color','k');
            fill([1:size(jTrqDatabt_m2,1) size(jTrqDatabt_m2,1):-1:1],...
                [jTrqDatabt_m2(:,ichan,idim,iset,iseq)+jTrqDatabt_std2(:,ichan,idim,iset,iseq) ;...
                jTrqDatabt_m2(end:-1:1,ichan,idim,iset,iseq)-jTrqDatabt_std2(end:-1:1,ichan,idim,iset,iseq)],...
                colb,'facealpha',0.4,'edgecolor','none');

            % PERT
            if any(iset == [1 2]) && (iseq > 3)
                for ipert = [1 5:8] % Don't plot inward perts
                    plot(jTrqDatapt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1.5);
        %             fill([1:size(jTrqDatapt_m2) size(jTrqDatapt_m2):-1:1],...
        %                 [jTrqDatapt_m2(:,ichan,idim,ipert,iset,iseq)+jTrqDatapt_std2(:,ichan,idim,ipert,iset,iseq) ;...
        %                 jTrqDatapt_m2(end:-1:1,ichan,idim,ipert,iset,iseq)-jTrqDatapt_std2(end:-1:1,ichan,idim,ipert,iset,iseq)],...
        %                 col(ipert,:),'facealpha',0.1,'edgecolor','none');
                end
            else
                for ipert = 1:8
                    plot(jTrqDatapt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1.5);
        %             fill([1:size(jTrqDatapt_m2) size(jTrqDatapt_m2):-1:1],...
        %                 [jTrqDatapt_m2(:,ichan,idim,ipert,iset,iseq)+jTrqDatapt_std2(:,ichan,idim,ipert,iset,iseq) ;...
        %                 jTrqDatapt_m2(end:-1:1,ichan,idim,ipert,iset,iseq)-jTrqDatapt_std2(end:-1:1,ichan,idim,ipert,iset,iseq)],...
        %                 col(ipert,:),'facealpha',0.1,'edgecolor','none');
                end
            end
            axis tight;
            
        end
    end

    % xy limits
    ymin = min(cellfun(@min,get(hplt,'ylim')));
    ymax = max(cellfun(@max,get(hplt,'ylim')));
    set(hplt,'xlim',[1 50],'ylim',[ymin ymax]);
    
    % xy labels
    set(hplt,'xtick',[]);
    xlabel(hplt(6),'TOR(PStart)-PEnd');
    xlabel(hplt(7),'PEnd-HSR');
    xlabel(hplt(8),'HSR-TOL');
    xlabel(hplt(9),'TOL-HSL');
    xlabel(hplt(10),'HSL-TOR');
    ylabel(hplt(1),chanXlbl{iplt});
    ylabel(hplt(6),chanXlbl{iplt});
    set(hplt([2:5 7:10]),'ytick',[],'ycolor',[1 1 1]);
    
    % Title
    suptitle(['Pert:' chanpTtl{iset} ' | Joint:' chanjTtl{iplt}]);

%     % Print figure
%     fullfilename = [pwd '\Figs\JTrq\' chanpTtl{iset} '_' chanjTtl{iplt}];
%     print(gcf,fullfilename,'-r300','-dpng')

end

end

%% Joint torques, on time averaged axis
%
% HipAP pos is (ante)flexion, neg is extension
% HipML pos is adduction, neg is abduction
% KneAP pos is extension, neg is flexion
% AnkAP pos is dorsiflexion, neg is plantarflexion
% AnkML pos is eversion, neg is inversion
% 
% Note that changes in hip joint angle can also occur due to pelvis rotation!

close all;
clear hplt;

sets = 4;%1:4;
for iset = sets

    
seqs = [1 2 4:6];
chandimL = [4 4 5 6 6 ; 1 2 1 1 2];
chandimR = [7 7 8 9 9 ; 1 2 1 1 2];
chanjTtl = {'HipAP','HipML','KneAP','AnkAP','AnkML'};
chanpTtl = {'SlowML','FastML','SlowAP','FastAP'};
chanXlbl = {{'Ext <> Flex','Nm/mgl'},{'Abd <> Add','Nm/mgl'},{'Flex <> Ext','Nm/mgl'},{'PFlex <> DFlex','Nm/mgl'},{'Inver <> Ever','Nm/mgl'}};

for iplt = 1:size(chandimL,2)
    
    figure('name',[chanpTtl{iset} '|' chanjTtl{iplt}],...
        'color',[1 1 1],...
        'position',[400 150 1080 720]);
    
    for iside = 0:1
        
        if iside == 0
            ichan = chandimL(1,iplt);
            idim = chandimL(2,iplt);
        else
            ichan = chandimR(1,iplt);
            idim = chandimR(2,iplt);
        end
        
        for iseq = seqs

            switch iseq
                case 1
                    hplt(1+5*iside) = subplot(2,6,1+6*iside); hold on;
                case 2
                    hplt(2+5*iside) = subplot(2,6,2+6*iside); hold on;
                case 4
                    hplt(3+5*iside) = subplot(2,6,3+6*iside); hold on;
                case 5
                    hplt(4+5*iside) = subplot(2,6,[4 5]+6*iside); hold on;
                case 6
                    hplt(5+5*iside) = subplot(2,6,6+6*iside); hold on;
            end

            % BASE
            plot(tDatabt_m2(:,iset,iseq),jTrqDatabt_m2(:,ichan,idim,iset,iseq),'-','linewidth',2,'color','k');
%             fill([1:size(jTrqDatabt_m2,1) size(jTrqDatabt_m2,1):-1:1],...
%                 [jTrqDatabt_m2(:,ichan,idim,iset,iseq)+jTrqDatabt_std2(:,ichan,idim,iset,iseq) ;...
%                 jTrqDatabt_m2(end:-1:1,ichan,idim,iset,iseq)-jTrqDatabt_std2(end:-1:1,ichan,idim,iset,iseq)],...
%                 colb,'facealpha',0.4,'edgecolor','none');

            % PERT
            if any(iset == [1 2]) && (iseq > 3)
                for ipert = [1 5:8] % Don't plot inward perts
                    plot(tDatapt_m2(:,ipert,iset,iseq),jTrqDatapt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1.5);
        %             fill([1:size(jTrqDatapt_m2) size(jTrqDatapt_m2):-1:1],...
        %                 [jTrqDatapt_m2(:,ichan,idim,ipert,iset,iseq)+jTrqDatapt_std2(:,ichan,idim,ipert,iset,iseq) ;...
        %                 jTrqDatapt_m2(end:-1:1,ichan,idim,ipert,iset,iseq)-jTrqDatapt_std2(end:-1:1,ichan,idim,ipert,iset,iseq)],...
        %                 col(ipert,:),'facealpha',0.1,'edgecolor','none');
                end
            else
                for ipert = 1:8
                    plot(tDatapt_m2(:,ipert,iset,iseq),jTrqDatapt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1.5);
%                     fill([1:size(jTrqDatapt_m2) size(jTrqDatapt_m2):-1:1],...
%                         [jTrqDatapt_m2(:,ichan,idim,ipert,iset,iseq)+jTrqDatapt_std2(:,ichan,idim,ipert,iset,iseq) ;...
%                         jTrqDatapt_m2(end:-1:1,ichan,idim,ipert,iset,iseq)-jTrqDatapt_std2(end:-1:1,ichan,idim,ipert,iset,iseq)],...
%                         col(ipert,:),'facealpha',0.1,'edgecolor','none');
                end
            end
            axis tight;
            
        end
    end

    % xy limits
    ymin = min(cellfun(@min,get(hplt,'ylim')));
    ymax = max(cellfun(@max,get(hplt,'ylim')));
%     set(hplt,'xlim',[1 50],'ylim',[ymin ymax]);
    set(hplt,'ylim',[ymin ymax]);
    
    % xy labels
    set(hplt,'xtick',[]);
    xlabel(hplt(6),'TOR(PStart)-PEnd');
    xlabel(hplt(7),'PEnd-HSR');
    xlabel(hplt(8),'HSR-TOL');
    xlabel(hplt(9),'TOL-HSL');
    xlabel(hplt(10),'HSL-TOR');
    ylabel(hplt(1),chanXlbl{iplt});
    ylabel(hplt(6),chanXlbl{iplt});
    set(hplt([2:5 7:10]),'ytick',[],'ycolor',[1 1 1]);
    
    % Title
    suptitle(['Pert:' chanpTtl{iset} ' | Joint:' chanjTtl{iplt}]);

    % Print figure
    fullfilename = [pwd '\Figs\JTrqTimeAx\' chanpTtl{iset} '_' chanjTtl{iplt}];
    print(gcf,fullfilename,'-r300','-dpng')

end

end


%% Joint power
%
% HipAP pos is (ante)flexion, neg is extension
% HipML pos is adduction, neg is abduction
% KneAP pos is extension, neg is flexion
% AnkAP pos is dorsiflexion, neg is plantarflexion
% AnkML pos is eversion, neg is inversion
% 
% Note that changes in hip joint angle can also occur due to pelvis rotation!

close all;
clear hplt;

sets = 1;%1:4;
for iset = sets
    
seqs = [1 2 4:6];
chandimL = [4 4 5 6 6 ; 1 2 1 1 2];
chandimR = [7 7 8 9 9 ; 1 2 1 1 2];
chanjTtl = {'HipAP','HipML','KneAP','AnkAP','AnkML'};
chanpTtl = {'SlowML','FastML','SlowAP','FastAP'};
chanXlbl = {'J s^-1 /(m*sqrt(g)*l^{1.5})','J s^-1/(m*sqrt(g)*l^{1.5})','J s^-1/(m*sqrt(g)*l^{1.5})','J/(m*sqrt(g)*l^{1.5})','J/(m*sqrt(g)*l^{1.5})'};

for iplt = 1:size(chandimL,2)
    
    figure('name',[chanpTtl{iset} '|' chanjTtl{iplt}],...
        'color',[1 1 1],...
        'position',[400 150 1080 720]);
    
    for iside = 0:1
        
        if iside == 0
            ichan = chandimL(1,iplt);
            idim = chandimL(2,iplt);
        else
            ichan = chandimR(1,iplt);
            idim = chandimR(2,iplt);
        end
        
        for iseq = seqs

            switch iseq
                case 1
                    hplt(1+5*iside) = subplot(2,6,1+6*iside); hold on;
                case 2
                    hplt(2+5*iside) = subplot(2,6,2+6*iside); hold on;
                case 4
                    hplt(3+5*iside) = subplot(2,6,3+6*iside); hold on;
                case 5
                    hplt(4+5*iside) = subplot(2,6,[4 5]+6*iside); hold on;
                case 6
                    hplt(5+5*iside) = subplot(2,6,6+6*iside); hold on;
            end

            % BASE
            plot(jPowDatabt_m2(:,ichan,idim,iset,iseq),'-','linewidth',2,'color','k');
            fill([1:size(jPowDatabt_m2,1) size(jPowDatabt_m2,1):-1:1],...
                [jPowDatabt_m2(:,ichan,idim,iset,iseq)+jPowDatabt_std2(:,ichan,idim,iset,iseq) ;...
                jPowDatabt_m2(end:-1:1,ichan,idim,iset,iseq)-jPowDatabt_std2(end:-1:1,ichan,idim,iset,iseq)],...
                colb,'facealpha',0.4,'edgecolor','none');

            % PERT
            if any(iset == [1 2]) && (iseq > 3)
                for ipert = [1 5:8]
                    plot(jPowDatapt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1.5);
        %             fill([1:size(jPowDatapt_m2) size(jPowDatapt_m2):-1:1],...
        %                 [jPowDatapt_m2(:,ichan,idim,ipert,iset,iseq)+jPowDatapt_std2(:,ichan,idim,ipert,iset,iseq) ;...
        %                 jPowDatapt_m2(end:-1:1,ichan,idim,ipert,iset,iseq)-jPowDatapt_std2(end:-1:1,ichan,idim,ipert,iset,iseq)],...
        %                 col(ipert,:),'facealpha',0.1,'edgecolor','none');
                end
            else
                for ipert = 1:8
                    plot(jPowDatapt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1.5);
        %             fill([1:size(jPowDatapt_m2) size(jPowDatapt_m2):-1:1],...
        %                 [jPowDatapt_m2(:,ichan,idim,ipert,iset,iseq)+jPowDatapt_std2(:,ichan,idim,ipert,iset,iseq) ;...
        %                 jPowDatapt_m2(end:-1:1,ichan,idim,ipert,iset,iseq)-jPowDatapt_std2(end:-1:1,ichan,idim,ipert,iset,iseq)],...
        %                 col(ipert,:),'facealpha',0.1,'edgecolor','none');
                end
            end
            
            axis tight;
            
        end
    end

    % xy limits
    ymin = min(cellfun(@min,get(hplt,'ylim')));
    ymax = max(cellfun(@max,get(hplt,'ylim')));
    set(hplt,'xlim',[1 50],'ylim',[ymin ymax]);
    
    % xy labels
    set(hplt,'xtick',[]);
    xlabel(hplt(6),'TOR(PStart)-PEnd');
    xlabel(hplt(7),'PEnd-HSR');
    xlabel(hplt(8),'HSR-TOL');
    xlabel(hplt(9),'TOL-HSL');
    xlabel(hplt(10),'HSL-TOR');
    ylabel(hplt(1),chanXlbl{iplt});
    ylabel(hplt(6),chanXlbl{iplt});
    set(hplt([2:5 7:10]),'ytick',[],'ycolor',[1 1 1]);
    
    % Title
    suptitle(['Pert:' chanpTtl{iset} ' | Joint:' chanjTtl{iplt}]);

%     % Print figure
%     fullfilename = [pwd '\Figs\JPow\' chanpTtl{iset} '_' chanjTtl{iplt}];
%     print(gcf,fullfilename,'-r300','-dpng')

end

end


%% Joint energy
%
% HipAP pos is (ante)flexion, neg is extension
% HipML pos is adduction, neg is abduction
% KneAP pos is extension, neg is flexion
% AnkAP pos is dorsiflexion, neg is plantarflexion
% AnkML pos is eversion, neg is inversion
% 
% Note that changes in hip joint angle can also occur due to pelvis rotation!

close all;
clear hplt;

sets = 1:4;
for iset = sets

    
seqs = [1 2 4:6];
chandimL = [4 4 5 6 6 ; 1 2 1 1 2];
chandimR = [7 7 8 9 9 ; 1 2 1 1 2];
chanjTtl = {'HipAP','HipML','KneAP','AnkAP','AnkML'};
chanpTtl = {'SlowML','FastML','SlowAP','FastAP'};
chanXlbl = {'J/mgl','J/mgl','J/mgl','J/mgl','J/mgl'};

for iplt = 1:size(chandimL,2)
    
    figure('name',[chanpTtl{iset} '|' chanjTtl{iplt}],...
        'color',[1 1 1],...
        'position',[400 150 1080 720]);
    
    for iside = 0:1
        
        if iside == 0
            ichan = chandimL(1,iplt);
            idim = chandimL(2,iplt);
        else
            ichan = chandimR(1,iplt);
            idim = chandimR(2,iplt);
        end
        
        for iseq = seqs

            switch iseq
                case 1
                    hplt(1+5*iside) = subplot(2,6,1+6*iside); hold on;
                case 2
                    hplt(2+5*iside) = subplot(2,6,2+6*iside); hold on;
                case 4
                    hplt(3+5*iside) = subplot(2,6,3+6*iside); hold on;
                case 5
                    hplt(4+5*iside) = subplot(2,6,[4 5]+6*iside); hold on;
                case 6
                    hplt(5+5*iside) = subplot(2,6,6+6*iside); hold on;
            end

            % BASE
%             plot(jEgyPosDatabt_m2(:,ichan,idim,iset,iseq),'-','linewidth',2,'color','k');
%             fill([1:size(jEgyPosDatabt_m2,1) size(jEgyPosDatabt_m2,1):-1:1],...
%                 [jEgyPosDatabt_m2(:,ichan,idim,iset,iseq)+jEgyPosDatabt_std2(:,ichan,idim,iset,iseq) ;...
%                 jEgyPosDatabt_m2(end:-1:1,ichan,idim,iset,iseq)-jEgyPosDatabt_std2(end:-1:1,ichan,idim,iset,iseq)],...
%                 colb,'facealpha',0.4,'edgecolor','none');
% 
%             plot(jEgyNegDatabt_m2(:,ichan,idim,iset,iseq),'-','linewidth',2,'color','k');
%             fill([1:size(jEgyNegDatabt_m2,1) size(jEgyNegDatabt_m2,1):-1:1],...
%                 [jEgyNegDatabt_m2(:,ichan,idim,iset,iseq)+jEgyNegDatabt_std2(:,ichan,idim,iset,iseq) ;...
%                 jEgyNegDatabt_m2(end:-1:1,ichan,idim,iset,iseq)-jEgyNegDatabt_std2(end:-1:1,ichan,idim,iset,iseq)],...
%                 colb,'facealpha',0.4,'edgecolor','none');
            
            plot(jEgyTotDatabt_m2(:,ichan,idim,iset,iseq),'-','linewidth',2,'color','k');
            fill([1:size(jEgyTotDatabt_m2,1) size(jEgyTotDatabt_m2,1):-1:1],...
                [jEgyTotDatabt_m2(:,ichan,idim,iset,iseq)+jEgyTotDatabt_std2(:,ichan,idim,iset,iseq) ;...
                jEgyTotDatabt_m2(end:-1:1,ichan,idim,iset,iseq)-jEgyTotDatabt_std2(end:-1:1,ichan,idim,iset,iseq)],...
                colb,'facealpha',0.4,'edgecolor','none');

            % PERT
%             if any(iset == [1 2])
%                 for ipert = [1 5:8]
% %                     plot(jEgyPosDatapt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1.5);
% %                     fill([1:size(jEgyPosDatapt_m2) size(jEgyPosDatapt_m2):-1:1],...
% %                         [jEgyPosDatapt_m2(:,ichan,idim,ipert,iset,iseq)+jEgyPosDatapt_std2(:,ichan,idim,ipert,iset,iseq) ;...
% %                         jEgyPosDatapt_m2(end:-1:1,ichan,idim,ipert,iset,iseq)-jEgyPosDatapt_std2(end:-1:1,ichan,idim,ipert,iset,iseq)],...
% %                         col(ipert,:),'facealpha',0.1,'edgecolor','none');
% 
% %                     plot(jEgyNegDatapt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1.5);
% %                     fill([1:size(jEgyNegDatapt_m2) size(jEgyNegDatapt_m2):-1:1],...
% %                         [jEgyNegDatapt_m2(:,ichan,idim,ipert,iset,iseq)+jEgyNegDatapt_std2(:,ichan,idim,ipert,iset,iseq) ;...
% %                         jEgyNegDatapt_m2(end:-1:1,ichan,idim,ipert,iset,iseq)-jEgyNegDatapt_std2(end:-1:1,ichan,idim,ipert,iset,iseq)],...
% %                         col(ipert,:),'facealpha',0.1,'edgecolor','none');
% 
%                     plot(jEgyTotDatapt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1.5);
% %                     fill([1:size(jEgyTotDatapt_m2) size(jEgyTotDatapt_m2):-1:1],...
% %                         [jEgyTotDatapt_m2(:,ichan,idim,ipert,iset,iseq)+jEgyTotDatapt_std2(:,ichan,idim,ipert,iset,iseq) ;...
% %                         jEgyTotDatapt_m2(end:-1:1,ichan,idim,ipert,iset,iseq)-jEgyTotDatapt_std2(end:-1:1,ichan,idim,ipert,iset,iseq)],...
% %                         col(ipert,:),'facealpha',0.1,'edgecolor','none');
%                 end
%             else
                for ipert = 1:8
%                     plot(jEgyPosDatapt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1.5);
%                     fill([1:size(jEgyPosDatapt_m2) size(jEgyPosDatapt_m2):-1:1],...
%                         [jEgyPosDatapt_m2(:,ichan,idim,ipert,iset,iseq)+jEgyPosDatapt_std2(:,ichan,idim,ipert,iset,iseq) ;...
%                         jEgyPosDatapt_m2(end:-1:1,ichan,idim,ipert,iset,iseq)-jEgyPosDatapt_std2(end:-1:1,ichan,idim,ipert,iset,iseq)],...
%                         col(ipert,:),'facealpha',0.1,'edgecolor','none');

%                     plot(jEgyNegDatapt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1.5);
%                     fill([1:size(jEgyNegDatapt_m2) size(jEgyNegDatapt_m2):-1:1],...
%                         [jEgyNegDatapt_m2(:,ichan,idim,ipert,iset,iseq)+jEgyNegDatapt_std2(:,ichan,idim,ipert,iset,iseq) ;...
%                         jEgyNegDatapt_m2(end:-1:1,ichan,idim,ipert,iset,iseq)-jEgyNegDatapt_std2(end:-1:1,ichan,idim,ipert,iset,iseq)],...
%                         col(ipert,:),'facealpha',0.1,'edgecolor','none');

                    plot(jEgyTotDatapt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1.5);
%                     fill([1:size(jEgyTotDatapt_m2) size(jEgyTotDatapt_m2):-1:1],...
%                         [jEgyTotDatapt_m2(:,ichan,idim,ipert,iset,iseq)+jEgyTotDatapt_std2(:,ichan,idim,ipert,iset,iseq) ;...
%                         jEgyTotDatapt_m2(end:-1:1,ichan,idim,ipert,iset,iseq)-jEgyTotDatapt_std2(end:-1:1,ichan,idim,ipert,iset,iseq)],...
%                         col(ipert,:),'facealpha',0.1,'edgecolor','none');
                end
%             end
            
            axis tight;
            
        end
    end

    % xy limits
    ymin = min(cellfun(@min,get(hplt,'ylim')));
    ymax = max(cellfun(@max,get(hplt,'ylim')));
    set(hplt,'xlim',[1 50],'ylim',[ymin ymax]);
    
    % xy labels
    set(hplt,'xtick',[]);
    xlabel(hplt(6),'TOR(PStart)-PEnd');
    xlabel(hplt(7),'PEnd-HSR');
    xlabel(hplt(8),'HSR-TOL');
    xlabel(hplt(9),'TOL-HSL');
    xlabel(hplt(10),'HSL-TOR');
    ylabel(hplt(1),chanXlbl{iplt});
    ylabel(hplt(6),chanXlbl{iplt});
    set(hplt([2:5 7:10]),'ytick',[],'ycolor',[1 1 1]);
    
    % Title
    suptitle(['Pert:' chanpTtl{iset} ' | Joint:' chanjTtl{iplt}]);

    % Print figure
    fullfilename = [pwd '\Figs\JEgy\' chanpTtl{iset} '_' chanjTtl{iplt}];
    print(gcf,fullfilename,'-r300','-dpng')

end

end


%% Joint angle vs joint torque
%
% HipAP pos is (ante)flexion, neg is extension
% HipML pos is adduction, neg is abduction
% KneAP pos is extension, neg is flexion
% AnkAP pos is dorsiflexion, neg is plantarflexion
% AnkML pos is eversion, neg is inversion
% 
% Note that changes in hip joint angle can also occur due to pelvis rotation!

close all;
clear hplt;

sets = 3;%1:4;
for iset = sets

    
seqs = [1 2 4:6];
chandimL = [4 4 5 6 6 ; 1 2 1 1 2];
chandimR = [7 7 8 9 9 ; 1 2 1 1 2];
chanjTtl = {'HipAP','HipML','KneAP','AnkAP','AnkML'};
chanpTtl = {'SlowML','FastML','SlowAP','FastAP'};
chanXlbl = {{'Ext <> Flex','rad'},{'Abd <> Add','rad'},{'Flex <> Ext','rad'},{'PFlex <> DFlex','rad'},{'Inver <> Ever','rad'}};
chanYlbl = {{'Ext <> Flex','Nm/mgl'},{'Abd <> Add','Nm/mgl'},{'Flex <> Ext','Nm/mgl'},{'PFlex <> DFlex','Nm/mgl'},{'Inver <> Ever','Nm/mgl'}};

for iplt = 1:size(chandimL,2)
    
    figure('name',[chanpTtl{iset} '|' chanjTtl{iplt}],...
        'color',[1 1 1],...
        'position',[400 150 1080 720]);
    
    for iside = 0:1
        
        if iside == 0
            ichan = chandimL(1,iplt);
            idim = chandimL(2,iplt);
        else
            ichan = chandimR(1,iplt);
            idim = chandimR(2,iplt);
        end
        
        for iseq = seqs

            hplt(1+iside) = subplot(1,2,1+iside); hold on;
            
            % BASE
            plot(jAngDatabt_m2(:,ichan,idim,iset,iseq),jTrqDatabt_m2(:,ichan,idim,iset,iseq),'-','linewidth',2,'color','k');
%             fill([1:size(jTrqDatabt_m2,1) size(jTrqDatabt_m2,1):-1:1],...
%                 [jTrqDatabt_m2(:,ichan,idim,iset,iseq)+jTrqDatabt_std2(:,ichan,idim,iset,iseq) ;...
%                 jTrqDatabt_m2(end:-1:1,ichan,idim,iset,iseq)-jTrqDatabt_std2(end:-1:1,ichan,idim,iset,iseq)],...
%                 colb,'facealpha',0.4,'edgecolor','none');

            % PERT
            if any(iset == [1 2])
                for ipert = [1 5:8]
                    plot(jAngDatapt_m2(:,ichan,idim,ipert,iset,iseq),jTrqDatapt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1.5);
        %             fill([1:size(jTrqDatapt_m2) size(jTrqDatapt_m2):-1:1],...
        %                 [jTrqDatapt_m2(:,ichan,idim,ipert,iset,iseq)+jTrqDatapt_std2(:,ichan,idim,ipert,iset,iseq) ;...
        %                 jTrqDatapt_m2(end:-1:1,ichan,idim,ipert,iset,iseq)-jTrqDatapt_std2(end:-1:1,ichan,idim,ipert,iset,iseq)],...
        %                 col(ipert,:),'facealpha',0.1,'edgecolor','none');
                end
            else
                for ipert = 1:8
                    plot(jAngDatapt_m2(:,ichan,idim,ipert,iset,iseq),jTrqDatapt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:),'linewidth',1.5);
        %             fill([1:size(jTrqDatapt_m2) size(jTrqDatapt_m2):-1:1],...
        %                 [jTrqDatapt_m2(:,ichan,idim,ipert,iset,iseq)+jTrqDatapt_std2(:,ichan,idim,ipert,iset,iseq) ;...
        %                 jTrqDatapt_m2(end:-1:1,ichan,idim,ipert,iset,iseq)-jTrqDatapt_std2(end:-1:1,ichan,idim,ipert,iset,iseq)],...
        %                 col(ipert,:),'facealpha',0.1,'edgecolor','none');
                end
            end
            axis tight;
            axis square;
        end
    end

    % xy limits
    xmin = min(cellfun(@min,get(hplt,'xlim')));
    xmax = max(cellfun(@max,get(hplt,'xlim')));
    ymin = min(cellfun(@min,get(hplt,'ylim')));
    ymax = max(cellfun(@max,get(hplt,'ylim')));
    set(hplt,'xlim',[xmin xmax],'ylim',[ymin ymax]);
    
    % xy labels
    xlabel(hplt(1),chanXlbl{iplt});
    ylabel(hplt(1),chanYlbl{iplt});
    xlabel(hplt(2),chanXlbl{iplt});
    ylabel(hplt(2),chanYlbl{iplt});
    
    % Title
    suptitle(['Pert:' chanpTtl{iset} ' | Joint:' chanjTtl{iplt}]);

%     % Print figure
%     fullfilename = [pwd '\Figs\JAvT\' chanpTtl{iset} '_' chanjTtl{iplt}];
%     print(gcf,fullfilename,'-r300','-dpng')

end

end



%% Pooled subjects, time series AP COM velocity, trunk cluster accel (for paper)
% (pinfoot paper)

close all;
% clear h1 h2 h3;
% 
% figure;
% 
% iset = 3;
% 
% for iplot = 1:4
% 
%     ncol = 7;
%     switch iplot
%         case 1
%             pltA = 1;
%             pltB = pltA + ncol;
% %             pltC = pltB + ncol;
%         case 2
%             pltA = [2 3];
%             pltB = pltA + ncol;
% %             pltC = pltB + ncol;
%         case 3
%             pltA = 4;
%             pltB = pltA + ncol;
% %             pltC = pltB + ncol;
%         case 4
%             pltA = [5 6 7];
%             pltB = pltA + ncol;
% %             pltC = pltB + ncol;
%     end
%     
%     % VELOCITY
%     h1(iplot) = subplot(2,ncol,pltA); hold on;
%     
%     % Baseline
%     fill([1:size(comDataDbt_m2,1) size(comDataDbt_m2,1):-1:1],...
%         [comDataDbt_m2(:,end,2,iset,iplot)+comDataDbt_std2(:,end,2,iset,iplot) ; comDataDbt_m2(end:-1:1,end,2,iset,iplot)-comDataDbt_std2(end:-1:1,end,2,iset,iplot)],...
%         colb,'facealpha',0.4,'edgecolor','none');
%     plot(comDataDbt_m2(:,end,2,iset,iplot),'color',colb,'linewidth',2);
% 
%     % Pert
%     for ipert = 1:8
% %         fill([1:size(comDataDpt_m2,1) size(comDataDpt_m2,1):-1:1],...
% %             [comDataDpt_m2(:,end,2,ipert,iplot)+comDataDpt_std2(:,end,2,ipert,iplot) ; comDataDpt_m2(end:-1:1,end,2,ipert,iplot)-comDataDpt_std2(end:-1:1,end,2,ipert,iplot)],...
% %             col(ipert,:),'facealpha',0.1,'edgecolor','none');
%         plot(comDataDpt_m2(:,end,2,ipert,iset,iplot),'color',col(ipert,:))
%     end
%     
%     % ANGLE
%     h2(iplot) = subplot(2,ncol,pltB); hold on;
%     
%     % Baseline
%     fill([1:length(cluDatabt_m2) length(cluDatabt_m2):-1:1],...
%     [cluDatabt_m2(:,8,1,iset,iplot)+cluDatabt_std2(:,8,1,iset,iplot) ; cluDatabt_m2(end:-1:1,8,1,iset,iplot)-cluDatabt_std2(end:-1:1,8,1,iset,iplot)],...
%     colb,'facealpha',0.4,'edgecolor','none');
%     plot(cluDatabt_m2(:,8,1,iset,iplot),'color',colb,'linewidth',2);
%     
% %     fill([1:length(cluDataDbt_m2) length(cluDataDbt_m2):-1:1],...
% %     [cluDataDbt_m2(:,8,1,iset,iplot)+cluDataDbt_std2(:,8,1,iset,iplot) ; cluDataDbt_m2(end:-1:1,8,1,iset,iplot)-cluDataDbt_std2(end:-1:1,8,1,iset,iplot)],...
% %     colb,'facealpha',0.4,'edgecolor','none');
% %     plot(cluDataDbt_m2(:,8,1,iset,iplot),'color',colb,'linewidth',2);
% 
% %     fill([1:length(cluDataDDbt_m2) length(cluDataDDbt_m2):-1:1],...
% %     [cluDataDDbt_m2(:,8,1,iset,iplot)+cluDataDDbt_std2(:,8,1,iset,iplot) ; cluDataDDbt_m2(end:-1:1,8,1,iset,iplot)-cluDataDDbt_std2(end:-1:1,8,1,iset,iplot)],...
% %     colb,'facealpha',0.4,'edgecolor','none');
% %     plot(cluDataDDbt_m2(:,8,1,iset,iplot),'color',colb,'linewidth',2);
%     
%     % Pert
%     for ipert = 1:8
%         fill([1:length(cluDatapt_m2) length(cluDatapt_m2):-1:1],...
%         [cluDatapt_m2(:,8,1,ipert,iset,iplot)+cluDatapt_std2(:,8,1,ipert,iset,iplot) ; cluDatapt_m2(end:-1:1,8,1,ipert,iset,iplot)-cluDatapt_std2(end:-1:1,8,1,ipert,iset,iplot)],...
%         col(ipert,:),'facealpha',0.1,'edgecolor','none');        
%         plot(cluDatapt_m2(:,8,1,ipert,iset,iplot),'color',col(ipert,:));
%         
% %         plot(cluDataDpt_m2(:,8,1,ipert,iset,iplot),'color',col(ipert,:));
%         
% %         fill([1:length(cluDataDDpt_m2) length(cluDataDDpt_m2):-1:1],...
% %         [cluDataDDpt_m2(:,8,1,ipert,iplot)+cluDataDDpt_std2(:,8,1,ipert,iplot) ; cluDataDDpt_m2(end:-1:1,8,1,ipert,iplot)-cluDataDDpt_std2(end:-1:1,8,1,ipert,iplot)],...
% %         col(ipert,:),'facealpha',0.1,'edgecolor','none');
% %         plot(cluDataDDpt_m2(:,8,1,ipert,iset,iplot),'color',col(ipert,:));
%     end
% 
% end
% 
% 
% set([h1 h2],'color','none');
% set([h1 h2],'xtick',[],'xticklabel',[],'xlim',[1 50],'xcolor','none');
% set([h1(2:end) h2(2:end)],'yticklabel',[],'ycolor','none');
% set([h1(1)],'ytick',-10:0.2:10);
% 
% set(h1,'ylim',[0.19 1.01]);
% set(h2,'ylim',[-0.151 0.151]); % Ang
% % set(h2,'ylim',[-1 1]); % Veloc
% % set(h2,'ylim',[-14 14]); % Accel
% 
% 
% ylabel(h1(1),{'AP COM velocity','[m/s]'});
% ylabel(h2(1),{'Trunk accel','[rad/s^2]'});
% 
% linkaxes(h1,'y');
% linkaxes(h2,'y');


% ####
% Alternative version, with 3 rows (COM vel, angle, angle acc)
% ####

clear h1 h2 h3
figure;

iset = 3;
seq2plot = [1 2 4 5 6];

for iplot = 1:5

    ncol = 7;
    switch iplot
        case 1
            pltA = 1;
            pltB = pltA + ncol;
            pltC = pltB + ncol;
        case 2
            pltA = [2 3];
            pltB = pltA + ncol;
            pltC = pltB + ncol;
        case 3
            pltA = 4;
            pltB = pltA + ncol;
            pltC = pltB + ncol;
        case 4
            pltA = [5 6];
            pltB = pltA + ncol;
            pltC = pltB + ncol;
        case 5
            pltA = 7;
            pltB = pltA + ncol;
            pltC = pltB + ncol;
    end
    
    % ###
    % VELOCITY
    % ###
    h1(iplot) = subplot(3,ncol,pltA); hold on;
    
    % Baseline
    fill([1:size(comDataDbt_m2,1) size(comDataDbt_m2,1):-1:1],...
        [comDataDbt_m2(:,end,2,iset,seq2plot(iplot))+comDataDbt_std2(:,end,2,iset,seq2plot(iplot)) ; comDataDbt_m2(end:-1:1,end,2,iset,seq2plot(iplot))-comDataDbt_std2(end:-1:1,end,2,iset,seq2plot(iplot))],...
        colb,'facealpha',0.4,'edgecolor','none');
    plot(comDataDbt_m2(:,end,2,iset,seq2plot(iplot)),'color',colb,'linewidth',2);

    % Pert
    for ipert = 1:8
%         fill([1:size(comDataDpt_m2,1) size(comDataDpt_m2,1):-1:1],...
%             [comDataDpt_m2(:,end,2,ipert,iplot)+comDataDpt_std2(:,end,2,ipert,iplot) ; comDataDpt_m2(end:-1:1,end,2,ipert,iplot)-comDataDpt_std2(end:-1:1,end,2,ipert,iplot)],...
%             col(ipert,:),'facealpha',0.1,'edgecolor','none');
        plot(comDataDpt_m2(:,end,2,ipert,iset,seq2plot(iplot)),'color',col(ipert,:))
    end
    
    % ###
    % ANGLE
    % ###
    h2(iplot) = subplot(3,ncol,pltB); hold on;
    
    % Baseline
    fill([1:length(cluDatabt_m2) length(cluDatabt_m2):-1:1],...
    [cluDatabt_m2(:,8,1,iset,seq2plot(iplot))+cluDatabt_std2(:,8,1,iset,seq2plot(iplot)) ; cluDatabt_m2(end:-1:1,8,1,iset,seq2plot(iplot))-cluDatabt_std2(end:-1:1,8,1,iset,seq2plot(iplot))],...
    colb,'facealpha',0.4,'edgecolor','none');
    plot(cluDatabt_m2(:,8,1,iset,seq2plot(iplot)),'color',colb,'linewidth',2);
    
%     fill([1:length(cluDataDbt_m2) length(cluDataDbt_m2):-1:1],...
%     [cluDataDbt_m2(:,8,1,iset,iplot)+cluDataDbt_std2(:,8,1,iset,iplot) ; cluDataDbt_m2(end:-1:1,8,1,iset,iplot)-cluDataDbt_std2(end:-1:1,8,1,iset,iplot)],...
%     colb,'facealpha',0.4,'edgecolor','none');
%     plot(cluDataDbt_m2(:,8,1,iset,iplot),'color',colb,'linewidth',2);

    
    % Pert
    for ipert = 1:8
%         fill([1:length(cluDatapt_m2) length(cluDatapt_m2):-1:1],...
%         [cluDatapt_m2(:,8,1,ipert,iset,iplot)+cluDatapt_std2(:,8,1,ipert,iset,iplot) ; cluDatapt_m2(end:-1:1,8,1,ipert,iset,iplot)-cluDatapt_std2(end:-1:1,8,1,ipert,iset,iplot)],...
%         col(ipert,:),'facealpha',0.1,'edgecolor','none');        
        plot(cluDatapt_m2(:,8,1,ipert,iset,seq2plot(iplot)),'color',col(ipert,:));
        
%         plot(cluDataDpt_m2(:,8,1,ipert,iset,iplot),'color',col(ipert,:));

    end

    
    % ###
    % ANGULAR VELOCITY
    % ###
    h3(iplot) = subplot(3,ncol,pltC); hold on;
    
    % Baseline
    fill([1:length(cluDataDbt_m2) length(cluDataDbt_m2):-1:1],...
    [cluDataDbt_m2(:,8,1,iset,seq2plot(iplot))+cluDataDbt_std2(:,8,1,iset,seq2plot(iplot)) ; cluDataDbt_m2(end:-1:1,8,1,iset,seq2plot(iplot))-cluDataDbt_std2(end:-1:1,8,1,iset,seq2plot(iplot))],...
    colb,'facealpha',0.4,'edgecolor','none');
    plot(cluDataDbt_m2(:,8,1,iset,seq2plot(iplot)),'color',colb,'linewidth',2);
    
    % Pert
    for ipert = 1:8
%         fill([1:length(cluDataDDpt_m2) length(cluDataDDpt_m2):-1:1],...
%         [cluDataDDpt_m2(:,8,1,ipert,iplot)+cluDataDDpt_std2(:,8,1,ipert,iplot) ; cluDataDDpt_m2(end:-1:1,8,1,ipert,iplot)-cluDataDDpt_std2(end:-1:1,8,1,ipert,iplot)],...
%         col(ipert,:),'facealpha',0.1,'edgecolor','none');
        plot(cluDataDpt_m2(:,8,1,ipert,iset,seq2plot(iplot)),'color',col(ipert,:));
    end
    
%     % ###
%     % ANGULAR ACCEL
%     % ###
%     h3(iplot) = subplot(3,ncol,pltC); hold on;
%     
%     % Baseline
%     fill([1:length(cluDataDDbt_m2) length(cluDataDDbt_m2):-1:1],...
%     [cluDataDDbt_m2(:,8,1,iset,seq2plot(iplot))+cluDataDDbt_std2(:,8,1,iset,seq2plot(iplot)) ; cluDataDDbt_m2(end:-1:1,8,1,iset,seq2plot(iplot))-cluDataDDbt_std2(end:-1:1,8,1,iset,seq2plot(iplot))],...
%     colb,'facealpha',0.4,'edgecolor','none');
%     plot(cluDataDDbt_m2(:,8,1,iset,seq2plot(iplot)),'color',colb,'linewidth',2);
%     
%     % Pert
%     for ipert = 1:8
% %         fill([1:length(cluDataDDpt_m2) length(cluDataDDpt_m2):-1:1],...
% %         [cluDataDDpt_m2(:,8,1,ipert,iplot)+cluDataDDpt_std2(:,8,1,ipert,iplot) ; cluDataDDpt_m2(end:-1:1,8,1,ipert,iplot)-cluDataDDpt_std2(end:-1:1,8,1,ipert,iplot)],...
% %         col(ipert,:),'facealpha',0.1,'edgecolor','none');
%         plot(cluDataDDpt_m2(:,8,1,ipert,iset,seq2plot(iplot)),'color',col(ipert,:));
%     end
%     
end

set([h1 h2 h3],'color','none');
set([h1 h2 h3],'xtick',[],'xticklabel',[],'xlim',[1 50],'xcolor','none');
set([h1(2:end) h2(2:end) h3(2:end)],'yticklabel',[],'ycolor','none');

set([h1(1)],'ytick',-10:0.2:10);

set(h1,'ylim',[0.19 1.01]);
set(h2,'ylim',[-0.151 0.151]); % Ang
% set(h3,'ylim',[-13.5 13.5]); % Accel
set(h3,'ylim',[-1 1]); % Veloc

ylabel(h1(1),{'AP COM velocity','[m/s]'});
ylabel(h2(1),{'Trunk ang. excursion','[rad]'});
ylabel(h3(1),{'Trunk ang. accel','[rad s^{-2}]'});

linkaxes(h1,'y');
linkaxes(h2,'y');
linkaxes(h3,'y');


%% Stacked gait phase durations 
% (pinfoot)

iset = 3;

tDatabe_m2c = zeros(size(tDatabe_m2(:,iset)));
tDatabe_m2c(3:end,:) = cumsum(tDatabe_m2(3:end,iset),1);
tDatape_m2c = zeros(size(tDatape_m2(:,:,iset)));
tDatape_m2c(3:end,:) = cumsum(tDatape_m2(3:end,:,iset),1);

% maxval = zeros(nevt,1);
% maxval(3:end) = max([tDatabe_m2(3:end,iset) tDatape_m2(3:end,:,iset)],[],2);
% maxcumval = cumsum( maxval );

% recwid = 0.5; % Rectangle width

figure; 
subplot(311);
hold on;
altord = [4:-1:1 6:9];
for ievt = 3:6
    
%     % Stacked
%     rectangle('position',[5,0+tDatabe_m2c(ievt-1),0.5,tDatabe_m2(ievt,iset)],...
%             'curvature',[0.5 0.5],'edgecolor',colb,'linewidth',2);
    
%     % Datapoints, y-axis is time axis
%     plot([5 5],[tDatabe_m2c(ievt) tDatabe_m2c(ievt)-tDatabe_std2(ievt,iset)],'-','color',colb,'linewidth',2); % 1-sided std
%     plot(5,tDatabe_m2c(ievt),'o','color',colb);
    
    % Datapoints, x-axis is time axis
    plot([tDatape_m2c(ievt,4:-1:1) tDatabe_m2c(ievt) tDatape_m2c(ievt,5:8)],1:9,'-k');
    plot([tDatabe_m2c(ievt)+tDatabe_std2(ievt,iset) tDatabe_m2c(ievt)-tDatabe_std2(ievt,iset)],[5 5],'-','color',colb,'linewidth',1); % 1-sided std
    plot(tDatabe_m2c(ievt),5,'o','color',colb,'markersize',4);
    
    for ipert = 1:8

%         % Stacked
%         rectangle('position',[altord(ipert),0+tDatape_m2c(ievt-1,ipert),0.5,tDatape_m2(ievt,ipert,iset)],...
%             'curvature',[0.5 0.5],'edgecolor',col(ipert,:),'linewidth',2);
        
        % Datapoints, y-axis is time axis
        plot([tDatape_m2c(ievt,ipert)+tDatape_std2(ievt,ipert,iset) tDatape_m2c(ievt,ipert)-tDatape_std2(ievt,ipert,iset)],[altord(ipert) altord(ipert)],'-','color',col(ipert,:),'linewidth',1);
        plot(tDatape_m2c(ievt,ipert),altord(ipert),lin2{ipert},'color',col(ipert,:),'markersize',4);
        
    end

end

% % For y-axis as time axis
% xlim([0.5 9.5]);
% ylim([0 1.5]);
% set(gca,'xcolor','none','color','none','ytick',0:0.3:1.5);
% ylabel('Cumulative duration')

% For x-axis as time axis
xlim([0 1.5]);
ylim([0.5 9.5]);
set(gca,'ycolor','none','color','none','xtick',0:0.3:1.5);
xlabel('Cumulative gait phase duration [s]')

%% Error with baseline (pinfoot)
% (area's)

close all;

figure;
altord = [4:-1:1 5:8];
offset = -0.05;

iTimeIdx = 1;

% For velocity
h1 = subplot(331); hold on;
plot(1:8,squeeze(comDataDDevpe_m2(iTimeIdx,end,2,altord,3)),'-k','linewidth',1.5);
for ipert = 1:8
%     plot(altord(ipert),comDataDDevpe_m2(1,end,2,ipert,3),lin2{ipert},'color',col(ipert,:));

%     errorbar(altord(ipert),comDataDDevpe_m2(2,end,2,ipert,3),comDataDDevpe_std2(2,end,2,ipert,3),'color',col(ipert,:));
    plot([altord(ipert) altord(ipert)]+offset,[comDataDDevpe_m2(iTimeIdx,end,2,ipert,3)+comDataDDevpe_std2(iTimeIdx,end,2,ipert,3) ...
        comDataDDevpe_m2(iTimeIdx,end,2,ipert,3)-comDataDDevpe_std2(iTimeIdx,end,2,ipert,3)],'color',col(ipert,:));
    plot(altord(ipert)+offset,comDataDDevpe_m2(iTimeIdx,end,2,ipert,3),lin2{ipert},'color',col(ipert,:),'linewidth',1,'markersize',4);
end


% For angle
h2 = subplot(334); hold on;
plot(1:8,squeeze(cluDataDevpe_m2(iTimeIdx,8,1,altord,3)),'-k','linewidth',1.5);
for ipert = 1:8
%     plot(altord(ipert),comDataDDevpe_m2(1,end,2,ipert,3),lin2{ipert},'color',col(ipert,:));

%     errorbar(altord(ipert),cluDataDevpe_m2(2,8,1,ipert,3),cluDataDevpe_std2(2,8,1,ipert,3),'color',col(ipert,:));
    plot([altord(ipert) altord(ipert)]+offset,[cluDataDevpe_m2(iTimeIdx,8,1,ipert,3)+cluDataDevpe_std2(iTimeIdx,8,1,ipert,3) ...
        cluDataDevpe_m2(iTimeIdx,8,1,ipert,3)-cluDataDevpe_std2(iTimeIdx,8,1,ipert,3)],'color',col(ipert,:));
    plot(altord(ipert)+offset,cluDataDevpe_m2(iTimeIdx,8,1,ipert,3),lin2{ipert},'color',col(ipert,:),'linewidth',1,'markersize',4);
end


% For angular accel
h2 = subplot(337); hold on;
plot(1:8,squeeze(cluDataDDDevpe_m2(iTimeIdx,8,1,altord,3)),'-k','linewidth',1.5);
for ipert = 1:8

%     errorbar(altord(ipert),cluDataDDDevpe_m2(2,8,1,ipert,3),cluDataDDDevpe_std2(2,8,1,ipert,3),'color',col(ipert,:));
    plot([altord(ipert) altord(ipert)]+offset,[cluDataDDDevpe_m2(iTimeIdx,8,1,ipert,3)+cluDataDDDevpe_std2(iTimeIdx,8,1,ipert,3) ...
        cluDataDDDevpe_m2(iTimeIdx,8,1,ipert,3)-cluDataDDDevpe_std2(iTimeIdx,8,1,ipert,3)],'color',col(ipert,:));
    plot(altord(ipert)+offset,cluDataDDDevpe_m2(iTimeIdx,8,1,ipert,3),lin2{ipert},'color',col(ipert,:),'linewidth',1,'markersize',4);
end

set([h1 h2],'xtick',[],'xcolor','none');

savefig(gcf,'Differences_NOPIN');

%% Grouped subjects: joint angles per gait phase
% 
% For ankle: negative is dorsiflexion...?

close all;
clear h;

iset = 1;

seqs = 3:6; % Sequences you want to plot
chans = [4 7];
% chans = [5 8];
% chans = [6 9];
dims = {1,1};

% chans = 4:9; % Channels you want to plot
% dims = {[1 2],[1],[1 2],[1 2],[1],[1 2]}; % Dimensions per channel

nseq = numel(seqs);
nchan = numel(chans);
ndim = sum(cellfun('prodofsize',dims));

figure;
for ichan = chans
    
    idxChan = find(chans == ichan);
    
    for idim = dims{idxChan}
        
        if idxChan > 1
            idxDim = find(dims{idxChan} == idim) + sum(cellfun('prodofsize',dims(1:idxChan-1)));
        else
            idxDim = find(dims{idxChan} == idim);
        end
        
        for iseq = seqs

            idxSeq = find(iseq == seqs);

            h( idxSeq+nseq*(idxDim-1) ) = subplot(ndim,nseq,idxSeq+nseq*(idxDim-1)); 
            hold on;

            % BASE
            plot(jAngDatabt_m2(:,ichan,idim,iset,iseq),'-','linewidth',2,'color',colb);
            fill([1:size(jAngDatabt_m2,1) size(jAngDatabt_m2,1):-1:1],...
                [jAngDatabt_m2(:,ichan,idim,iset,iseq)+jAngDatabt_std2(:,ichan,idim,iset,iseq) ;...
                jAngDatabt_m2(end:-1:1,ichan,idim,iset,iseq)-jAngDatabt_std2(end:-1:1,ichan,idim,iset,iseq)],...
                colb,'facealpha',0.4,'edgecolor','none');

            % PERT
            for ipert = 1:8
                plot(jAngDatapt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:));
%                 fill([1:size(jAngDatapt_m2) size(jAngDatapt_m2):-1:1],...
%                     [jAngDatapt_m2(:,ichan,idim,ipert,iset,iseq)+jAngDatapt_std2(:,ichan,idim,ipert,iset,iseq) ;...
%                     jAngDatapt_m2(end:-1:1,ichan,idim,ipert,iset,iseq)-jAngDatapt_std2(end:-1:1,ichan,idim,ipert,iset,iseq)],...
%                     col(ipert,:),'facealpha',0.1,'edgecolor','none');
            end
            
        end
        
    end
end

set(h,'ylim',[-0.5 0.5]);
linkaxes(h,'xy');


%% Grouped subjects: joint torques per gait phase

close all;
clear h;

iset = 3;

seqs = 3:6; % Sequences you want to plot
chans = [4 7];
% chans = [5 8];
% chans = [6 9];
dims = {1,1};

% chans = 4:9; % Channels you want to plot
% dims = {[1 2],[1],[1 2],[1 2],[1],[1 2]}; % Dimensions per channel

nseq = numel(seqs);
nchan = numel(chans);
ndim = sum(cellfun('prodofsize',dims));

figure;
for ichan = chans
    
    idxChan = find(chans == ichan);
    
    for idim = dims{idxChan}
        
        if idxChan > 1
            idxDim = find(dims{idxChan} == idim) + sum(cellfun('prodofsize',dims(1:idxChan-1)));
        else
            idxDim = find(dims{idxChan} == idim);
        end
        
        for iseq = seqs

            idxSeq = find(iseq == seqs);

            h( idxSeq+nseq*(idxDim-1) ) = subplot(ndim,nseq,idxSeq+nseq*(idxDim-1)); 
            hold on;

            % BASE
            plot(jTrqDatabt_m2(:,ichan,idim,iset,iseq),'-','linewidth',2,'color',colb);
            fill([1:size(jTrqDatabt_m2,1) size(jTrqDatabt_m2,1):-1:1],...
                [jTrqDatabt_m2(:,ichan,idim,iset,iseq)+jTrqDatabt_std2(:,ichan,idim,iset,iseq) ;...
                jTrqDatabt_m2(end:-1:1,ichan,idim,iset,iseq)-jTrqDatabt_std2(end:-1:1,ichan,idim,iset,iseq)],...
                colb,'facealpha',0.4,'edgecolor','none');

            % PERT
            for ipert = 1:8
                plot(jTrqDatapt_m2(:,ichan,idim,ipert,iset,iseq),'-','color',col(ipert,:));
%                 fill([1:size(jTrqDatapt_m2) size(jTrqDatapt_m2):-1:1],...
%                     [jTrqDatapt_m2(:,ichan,idim,ipert,iset,iseq)+jTrqDatapt_std2(:,ichan,idim,ipert,iset,iseq) ;...
%                     jTrqDatapt_m2(end:-1:1,ichan,idim,ipert,iset,iseq)-jTrqDatapt_std2(end:-1:1,ichan,idim,ipert,iset,iseq)],...
%                     col(ipert,:),'facealpha',0.1,'edgecolor','none');
            end
            
        end
        
    end
end

set(h,'xticklabel',[],'xlim',[0 nsmpl],'ylim',[-0.1 0.1]);

firstCol = 1:nseq:ndim*nseq;
otherCol = setdiff(1:ndim*nseq,firstCol);
lastRow = 1+(ndim-1)*nseq:ndim*nseq;
otherRow = 1:(ndim-1)*nseq;

for iplot = firstCol
    ylabel(h(iplot),'Torque [a.u.]');
end
set(h(otherCol),'yticklabel',[]);

for iplot = lastRow
    seq = find(lastRow==iplot);
    switch seq
        case 1
            xlabel(h(lastRow(seq)),'TOR-HSR');
        case 2
            xlabel(h(lastRow(seq)),'HSR-TOL');
        case 3
            xlabel(h(lastRow(seq)),'TOL-HSL');
        case 4
            xlabel(h(lastRow(seq)),'HSL-TOR');
    end
    
end

linkaxes(h,'xy');




%% Individual subjects: joint angles per gait phase
close all

iset = 4;
% isubj = 1;

for isubj = 1:10
    
seqs = 3:6; % Sequences you want to plot
chans = [4 7];
% chans = [5 8];
% chans = [6 9];
dims = {1,1};

% chans = 4:9; % Channels you want to plot
% dims = {[1 2],[1],[1 2],[1 2],[1],[1 2]}; % Dimensions per channel

nseq = numel(seqs);
nchan = numel(chans);
ndim = sum(cellfun('prodofsize',dims));

figure;
for ichan = chans
    
    idxChan = find(chans == ichan);
    
    for idim = dims{idxChan}
        
        if idxChan > 1
            idxDim = find(dims{idxChan} == idim) + sum(cellfun('prodofsize',dims(1:idxChan-1)));
        else
            idxDim = find(dims{idxChan} == idim);
        end
        
        for iseq = seqs

            idxSeq = find(iseq == seqs);

            h( idxSeq+nseq*(idxDim-1) ) = subplot(ndim,nseq,idxSeq+nseq*(idxDim-1)); 
            hold on;

            % BASE
            plot(jAngDatabt_m1(:,ichan,idim,iset,isubj,iseq),'-','linewidth',2,'color',colb);
            fill([1:size(jAngDatabt_m1,1) size(jAngDatabt_m1,1):-1:1],...
                [jAngDatabt_m1(:,ichan,idim,iset,isubj,iseq)+jAngDatabt_std1(:,ichan,idim,iset,isubj,iseq) ;...
                jAngDatabt_m1(end:-1:1,ichan,idim,iset,isubj,iseq)-jAngDatabt_std1(end:-1:1,ichan,idim,iset,isubj,iseq)],...
                colb,'facealpha',0.4,'edgecolor','none');

            % PERT
            for ipert = 1:8
                plot(jAngDatapt_m1(:,ichan,idim,ipert,iset,isubj,iseq),'-','color',col(ipert,:));
                fill([1:size(jAngDatapt_m1) size(jAngDatapt_m1):-1:1],...
                    [jAngDatapt_m1(:,ichan,idim,ipert,iset,isubj,iseq)+jAngDatapt_std1(:,ichan,idim,ipert,iset,isubj,iseq) ;...
                    jAngDatapt_m1(end:-1:1,ichan,idim,ipert,iset,isubj,iseq)-jAngDatapt_std1(end:-1:1,ichan,idim,ipert,iset,isubj,iseq)],...
                    col(ipert,:),'facealpha',0.1,'edgecolor','none');
            end
            
        end
        
    end
end

set(h,'ylim',[-1 1]);
linkaxes(h,'xy');


end



%% For Herman, WeRob 2016, 1 subject, COM foot vs COM velocity

close all

isubj = 4;
iseq = 3;

figure; 
h1=subplot(121); hold on;

plot([-10 10],[-10 10].*w0inv(isubj),'-','color',[1 0.5 1],'linewidth',1.5);
plot(comDataDbt_m1(1,end,2,3,isubj,iseq),comDatabt_m1(1,7,2,3,isubj,iseq),'o','color',colb,'linewidth',1.5);
% plot(comDataDbt_m1(1,end,2,3,isubj,iseq),comDatabt_m1(1,4,2,3,isubj,iseq),'o','color',colb);

for ipert = 1:8
    
    plot(comDataDpt_m1(1,end,2,ipert,3,isubj,iseq),comDatapt_m1(1,7,2,ipert,3,isubj,iseq),lin2{ipert},'color',col(ipert,:),'linewidth',1.5);
%     plot(comDataDpt_m1(1,end,2,ipert,3,isubj,iseq),comDatapt_m1(1,4,2,ipert,3,isubj,iseq),lin2{ipert},'color',col(ipert,:));
    
end

set(h1,'xlim',[0.2 1],'ylim',[-0.1 0.3]);
axis square;

% % Averages
% close all
% 
% iseq = 3;
% 
% figure; 
% h1=subplot(121); hold on;
% 
% plot([-10 10],[-10 10].*mean(w0inv(:)),'-','color',[1 0.5 1],'linewidth',1.5);
% plot(comDataDbt_m2(1,end,2,3,iseq),comDatabt_m2(1,7,2,3,iseq),'o','color',colb,'linewidth',1.5);
% % plot(comDataDbt_m1(1,end,2,3,isubj,iseq),comDatabt_m1(1,4,2,3,isubj,iseq),'o','color',colb);
% 
% for ipert = 1:8
%     
%     plot(comDataDpt_m2(1,end,2,ipert,3,iseq),comDatapt_m2(1,7,2,ipert,3,iseq),lin2{ipert},'color',col(ipert,:),'linewidth',1.5);
% %     plot(comDataDpt_m1(1,end,2,ipert,3,isubj,iseq),comDatapt_m1(1,4,2,ipert,3,isubj,iseq),lin2{ipert},'color',col(ipert,:));
%     
% end
% 
% set(h1,'xlim',[0.2 1],'ylim',[-0.1 0.3]);
% axis square;


%%
% close all
% 
% iseq = 3;
% 
% figure; 
% subplot(121); hold on;
% 
% plot(comDatabt_m2(1,7,1,3,iseq),comDatabt_m2(1,7,2,3,iseq),'o','color',colb);
% plot(comDatabt_m2(1,4,1,3,iseq),comDatabt_m2(1,4,2,3,iseq),'o','color',colb);
% 
% for ipert = 1:8
%     
%     plot(comDatapt_m2(1,7,1,ipert,3,iseq),comDatapt_m2(1,7,2,ipert,3,iseq),lin2{ipert},'color',col(ipert,:));
%     plot(comDatapt_m2(1,4,1,ipert,3,iseq),comDatapt_m2(1,4,2,ipert,3,iseq),lin2{ipert},'color',col(ipert,:));
%     
% end





