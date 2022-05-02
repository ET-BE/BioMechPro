function [Datastr] = getmotioncycle(Datastr,angleDatafrom)
%% Motion Cycle detection based on joint angle data
% INPUT)
% Datastr : str
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
%
% phasevectorlabel : cell, corresponding with the columns in phasevector
% Contains the names of the joints
%
% phaseidx : structure with the indices of the start and end of each cycle
%
% NOTES)

% The function assumes that the data is smooth and properly filtered !
% Motion cycle is defined in the sagittal plane, i.e. Flexion extension (FE)
%
%
% TODO

% Turn off warning signal:findpeaks:largeMinPeakHeight. Some motion tasks,
% e.g. gait, wont have the min. height
warning('off','signal:findpeaks:largeMinPeakHeight')
%% Check on input data and number of outputs
if strcmp(angleDatafrom,'Marker')
    try
        % nJoints=size(Datastr.Marker.JointAngData,2);
        ch_jnt=[8,9]; %Only knee and ankle right
        nJoints=ch_jnt;
        anglelabel=Datastr.Marker.JointAngDataLabel;
        angleData=Datastr.Marker.JointAngData(:,:,1)./pi*180; %Only in FE angle
        fs=Datastr.Marker.MarkerFrameRate;
    catch
        error('getmotioncycle:angledata','No field Marker.JointAngData found in input data structure');
    end
    
    
    %% Get
    
    %1 seg distance
    for iCh=nJoints
        
        clear p_max p_min samp_max samp_min
        %     [p_max,samp_max]=findpeaks(angleData(:,iCh));
        peakdist=1; %1 seg distance Assume movement to be slower than 1Hz
        switch anglelabel{iCh}
            case 'KJCR'
                 %Gait (0 to-50deg)
                %Calf contrac. (0 to -10deg)
                %Squats (0 to -110deg) min -80deg -> of interest for squats
                %interesting for squats
                %                 peakheight=-(-80); %as -angleData in findpeak
                %                 [p_min,samp_min]=findpeaks(-angleData(:,iCh),'MinPeakDistance',peakdist*fs,'MinPeakHeight',peakheight); % 1s distance 1s*fs (100 if fs=100)
                %                  if ~isempty(samp_min) && size(samp_min,1)>2
                %                     Kneestart=samp_min(1:end-1);
                %                     Kneeend=samp_min(2:end)-1;
                %                 else
                %                     Kneestart=[];
                %                     Kneeend=[];
                %                 end
                peakheight=-20;
                [p_max,samp_max]=findpeaks(angleData(:,iCh),'MinPeakDistance',peakdist*fs,'MinPeakHeight',peakheight); % 1s distance 1s*fs (100 if fs=100)
                p_min=[]; samp_min=[];
                if ~isempty(samp_max) && size(samp_max,1)>2
                    Kneestart=samp_max(1:end-1);
                    Kneeend=samp_max(2:end)-1;
                else
                    Kneestart=[];
                    Kneeend=[];
                end
            case 'AJCR'
                %Gait (25 to -20deg)
                %Calf contrac. (10 to -40deg) min -20-> of interest for
                %calf contrc
                %Squats (0 to -110deg) min -80deg
                %interesting for squats
                %                 peakheight=-(-20); %as -angleData in findpeak
                %                 [p_min,samp_min]=findpeaks(-angleData(:,iCh),'MinPeakDistance',peakdist*fs,'MinPeakHeight',peakheight); % 1s distance 1s*fs (100 if fs=100)
                %                 if ~isempty(samp_min) && size(samp_min,1)>2
                %                     Anklestart=samp_min(1:end-1);
                %                     Ankleend=samp_min(2:end)-1;
                %                 else
                %                     Anklestart=[];
                %                     Ankleend=[];
                %                 end
                peakheight=-(-15); %as -angleData in findpeak
                [p_min,samp_min]=findpeaks(-angleData(:,iCh),'MinPeakDistance',peakdist*fs,'MinPeakHeight',peakheight); % 1s distance 1s*fs (100 if fs=100)
                p_max=[];  samp_max=[];
                
                if ~isempty(samp_min) && size(samp_min,1)>=2                    
                    samp_min_aux=round((samp_min(1:end-1)+ (samp_min(2:end)-1))/2);
                    plot(n,angleData(:,iCh),samp_min_aux,angleData(samp_min_aux,iCh),'v');
                    samp_min=samp_min_aux;
                    p_min=-angleData(samp_min,iCh);
                    
                    Anklestart=samp_min(1:end-1);
                    Ankleend=samp_min(2:end)-1;
                else
                    Anklestart=[];
                    Ankleend=[];
                end
                %otherwise include more joints
        end
        
        p_min=-p_min;
        %Plot with peaks
        figure('name', ['Trial ' num2str(Datastr.Info.Trial) '| ' anglelabel{iCh}],'NumberTitle','off')
        n=1:length(angleData(:,iCh));
        %     x=plot(n,angleData(:,iCh),samp_max,p_max,'^',samp_min,p_min,'v');
        plot(n,angleData(:,iCh),samp_min,p_min,'v');
    end
    % plot(n,angleData(:,iCh),Kneestart,p_min(1:end-1),'r*',Kneeend,p_min(2:end),'g*');
    
    
    
elseif strcmp(angleDatafrom,'IMU')
    try
        % nJoints=size(Datastr.Marker.JointAngData,2);
        ch_jnt=[10,13]; %Only knee and ankle right
        nJoints=ch_jnt;
        anglelabel=Datastr.IMU.IMUDataLabel;
        angleData=Datastr.IMU.IMUData(:,:); %Only in FE angle
        fs=Datastr.IMU.IMUFrameRate;
    catch
        error('getmotioncycle:angledata','No field IMU.IMUData found in input data structure');
    end
    
    %% Get
    %Based on joint angles
    %1 seg distance
    for iCh=nJoints
        
        clear p_max p_min samp_max samp_min
        %     [p_max,samp_max]=findpeaks(angleData(:,iCh));
        peakdist=1; %1 seg distance Assume movement to be slower than 1Hz
        switch anglelabel{iCh}
            case 'knee_angle_r'
                %Gait (0 to-50deg)
                %Calf contrac. (0 to -10deg)
                %Squats (0 to -110deg) min -80deg -> of interest for squats
                %interesting for squats
                %                 peakheight=-(-80); %as -angleData in findpeak
                %                 [p_min,samp_min]=findpeaks(-angleData(:,iCh),'MinPeakDistance',peakdist*fs,'MinPeakHeight',peakheight); % 1s distance 1s*fs (100 if fs=100)
                %                  if ~isempty(samp_min) && size(samp_min,1)>2
                %                     Kneestart=samp_min(1:end-1);
                %                     Kneeend=samp_min(2:end)-1;
                %                 else
                %                     Kneestart=[];
                %                     Kneeend=[];
                %                 end
                peakheight=-20;
                [p_max,samp_max]=findpeaks(angleData(:,iCh),'MinPeakDistance',peakdist*fs,'MinPeakHeight',peakheight); % 1s distance 1s*fs (100 if fs=100)
                p_min=[]; samp_min=[];
                if ~isempty(samp_max) && size(samp_max,1)>2
                    Kneestart=samp_max(1:end-1);
                    Kneeend=samp_max(2:end)-1;
                else
                    Kneestart=[];
                    Kneeend=[];
                end
            case 'ankle_angle_r'
                %Gait (25 to -20deg)
                %Calf contrac. (10 to -40deg) min -20-> of interest for
                %calf contrc
                %Squats (0 to -110deg) min -80deg
                %interesting for squats
                %                 peakheight=-(-20); %as -angleData in findpeak
                %                 [p_min,samp_min]=findpeaks(-angleData(:,iCh),'MinPeakDistance',peakdist*fs,'MinPeakHeight',peakheight); % 1s distance 1s*fs (100 if fs=100)
                %                 if ~isempty(samp_min) && size(samp_min,1)>2
                %                     Anklestart=samp_min(1:end-1);
                %                     Ankleend=samp_min(2:end)-1;
                %                 else
                %                     Anklestart=[];
                %                     Ankleend=[];
                %                 end
                peakheight=-(-15); %as -angleData in findpeak
                [p_min,samp_min]=findpeaks(-angleData(:,iCh),'MinPeakDistance',peakdist*fs,'MinPeakHeight',peakheight); % 1s distance 1s*fs (100 if fs=100)
                p_max=[];  samp_max=[];
                
                if ~isempty(samp_min) && size(samp_min,1)>=2                    
                    samp_min_aux=round((samp_min(1:end-1)+ (samp_min(2:end)-1))/2);
%                     n=1:length(angleData(:,iCh));
%                     plot(n,angleData(:,iCh),samp_min_aux,angleData(samp_min_aux,iCh),'v');
                    samp_min=samp_min_aux;
                    p_min=-angleData(samp_min,iCh);
                    
                    Anklestart=samp_min(1:end-1);
                    Ankleend=samp_min(2:end)-1;
                else
                    Anklestart=[];
                    Ankleend=[];
                end
                %otherwise include more joints
        end
        
        p_min=-p_min;
        %Plot with peaks
        figure('name', ['Trial ' num2str(Datastr.Info.Trial) '| ' anglelabel{iCh}],'NumberTitle','off')
        n=1:length(angleData(:,iCh));
        plot(n,angleData(:,iCh),samp_max,p_max,'^',samp_min,p_min,'v');
    end
    % plot(n,angleData(:,iCh),Kneestart,p_min(1:end-1),'r*',Kneeend,p_min(2:end),'g*');
    
    %Based on contact data from imu
    
    if isfield(Datastr.IMU,'IMUContactRFoot')
        peakdist=0.5*fs; %1 seg distance Assume movement to be slower than 1Hz
        %         [p_min,samp_min]=findpeaks(-Datastr.IMU.IMUContactRFoot); % 1s distance 1s*fs (100 if fs=100)
        samp_min_start=find(diff(Datastr.IMU.IMUContactRFoot)==1);
        samp_min_end=find(diff(Datastr.IMU.IMUContactRFoot)==-1);
        %Check order of starting
        f= figure('name', ['Trial ' num2str(Datastr.Info.Trial) '| ' 'ContactData'],'NumberTitle','off');
        n=1:length(Datastr.IMU.IMUContactRFoot);
        plot(n./fs,Datastr.IMU.IMUContactRFoot,samp_min_start./fs,ones(1,length(samp_min_start)),'*',samp_min_end./fs,ones(1,length(samp_min_end)),'o'); hold on
        t=[1:length(Datastr.EMG.normEMGData(:,1))]./(Datastr.EMG.EMGFrameRate);
        plot(t,Datastr.EMG.normEMGData(:,1)*4)
        
        prompt = 'Reverse start and end in contact data?';
        title = 'Input';
        dims = [1 35];
        definput = {'0','hsv'};
        answer = inputdlg(prompt,title,dims,definput);
        x=str2num(answer{1});
        
        if x==1
            samp_min_end=find(diff(Datastr.IMU.IMUContactRFoot)==1);
            samp_min_start=find(diff(Datastr.IMU.IMUContactRFoot)==-1);
            close(f)
                     figure('name', ['Trial ' num2str(Datastr.Info.Trial) '| ' 'ContactData'],'NumberTitle','off')
        n=1:length(Datastr.IMU.IMUContactRFoot);
        b=abs(Datastr.IMU.IMUContactRFoot-1);
        plot(n./fs,b,samp_min_start./fs,ones(1,length(samp_min_start)),'*',samp_min_end./fs,ones(1,length(samp_min_end)),'o'); hold on
        t=[1:length(Datastr.EMG.normEMGData(:,1))]./(Datastr.EMG.EMGFrameRate);
        plot(t,Datastr.EMG.normEMGData(:,1)*4)
        end
        
        if ~isempty(samp_min_start) && ~isempty(samp_min_end)
            %Ensure equal number of rising and falling edges
            N=min(length(samp_min_start),length(samp_min_end));
            samp_min_start=samp_min_start(1:N);
            samp_min_end=samp_min_end(1:N);
            %Ensure we start with a rising edge
            if samp_min_start(1)>samp_min_end(1)
                samp_min_start=samp_min_start(1:end-1);
                samp_min_end=samp_min_end(2:end);
            end
            
            valid=(samp_min_end-samp_min_start)>=peakdist;
            Stanceend=samp_min_end(valid);
            Stancestart=samp_min_start(valid);
            samp_min=[Stancestart  Stanceend];
            p_min=Datastr.IMU.IMUContactRFoot(samp_min);
            p_max=[];  samp_max=[];
        else
            Stancestart=[];
            Stanceend=[];
        end
        %Debugging
        %Debugging
%               figure('name', ['Trial ' num2str(Datastr.Info.Trial) '| ' 'ContactData'],'NumberTitle','off')
%         n=1:length(Datastr.IMU.IMUContactRFoot);
%         plot(n./fs,Datastr.IMU.IMUContactRFoot,samp_min_start./fs,ones(1,length(samp_min_start))*100,'*',samp_min_end./fs,ones(1,length(samp_min_end))*100,'o'); hold on
%         plot(n./fs,Datastr.IMU.IMUContactRFoot,samp_min(:,1)./fs,ones(1,length(samp_min(:,1))),'*',samp_min(:,2)./fs,ones(1,length(samp_min(:,2))),'o'); hold on
%         t=[1:length(Datastr.Force.ForceData(:,9))]./(Datastr.Force.ForceFrameRate);
%         plot(t,Datastr.Force.ForceData(:,9))
        %Plot with peaks
%         figure('name', ['Trial ' num2str(Datastr.Info.Trial) '| ' 'ContactData'],'NumberTitle','off')
%         n=1:length(Datastr.IMU.IMUContactRFoot);
%         plot(n./fs,Datastr.IMU.IMUContactRFoot,samp_max./fs,p_max,'^',samp_min./fs,p_min,'v'); hold on
%         n=1:length(Datastr.IMU.IMUData(:,nJoints));
%         plot(n/fs,Datastr.IMU.IMUData(:,nJoints(1)),n/fs,Datastr.IMU.IMUData(:,nJoints(2)));
%         phaseidx.IMUR = [Stancestart Stanceend];
%     end
    
    
end
% Turn on warning signal:findpeaks:largeMinPeakHeight.
warning('on','signal:findpeaks:largeMinPeakHeight')

phasevectorlabel = {'KneeRcycle','AnkleRcycle','IMURstance'};
phaseidx.KneeR = [Kneestart Kneeend];
phaseidx.AnkleR = [Anklestart Ankleend];
phaseidx.IMUR = [Stancestart  Stanceend]
%% Generate output
Datastr.Event.MotionCycleLabel = phasevectorlabel;
Datastr.Event.MotionCycleIdx = phaseidx;
end
