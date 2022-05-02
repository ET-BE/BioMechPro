function [Datastr] = E2_procEMG(Datastr,EMGs_Input,setOutliersTo1)
% gBMPDynUI EMGs_Input=1; setOutliersTo1=1;
%EMGs_Input={'TAL','GML','RFL'};
% Detrend, rectify, and low pass filter EMG
%
% INPUT)
% Datastr: structure, with at least the fields
%     .EMG.EMGData
%     .EMG.EMGFrameRate
% setOutliersTo1= %=1; if EMG surpasses the maximum value found in MVC, it is set to 1
%                 %=0; values surpassing the max emg value will have values higher than 1
%
% OUTPUT)
% Datastr: structure with no new added fields
%     .EMG.filtEMGData
%     .EMG.normEMGData
%     .EMG.filtEMGtime
%
% NOTES)
% Function based on MOtion data elaboration TOolbox for
% NeuroMusculoSkeletal applications (MOtoNMS).
% Copyright (C) 2012-2014 Alice Mantoan, Monica Reggiani

% TODO
%This module could be more efficient, loading the trial to be normalized in
%each trial processing from BioMechPro. Nevbertheless, to keep MOtoNMS
%functions the same, all the processing is done when we loading the first
%trial in BioMechpro. All trials are processed in the first loading. The next loadings are skipped 
%(as they were already proccessed during the first trial. Read the code!)


%% Check input

try %load flagEMG.
    cd(Datastr.Info.SubjRoot)
    load('flagEMG.mat')
catch %flagEMG won't exist for the first trial of a subject, create a flag
    flagEMG=0;
    save('flagEMG.mat','flagEMG')
end

if flagEMG %If flagEMG==1, all trials for this subject were already filtered and normalized. We dont need to proccess the subject EMGs anymore
    if isfield(Datastr,'EMG')
        fprintf(['Trial ' Datastr.Info.Trial ' processed (flagEMG=1).\n ']) % All trials were processsed during the first trial        
    else
        warning(['No EMG data in file ' Datastr.Info.Trial '. Skipping.'])
    end
    Datastr=[];
    return
else
    
    %% Do Stuff
    %Input BioMechPro
    
    % %% ------------------------------------------------------------------------
    % %                      START/STOP COMPUTATION
    % %--------------------------------------------------------------------------

    maxemgPath=[Datastr.Info.SubjRoot filesep 'EMG' filesep 'maxemg'];
    mkdir(maxemgPath);
    listing = dir('*.mat');
    eraseIdx=[];
    iDx=[];
    fprintf('\n#### Loading trials to be normalized... ####\n')
    for i=1:numel(listing) %List .mat files in folder of subject
        matList{i}=listing(i).name;
        if ~strcmp(matList{i},'flagEMG.mat')
        load(matList{i});
        if isfield(Datastr,'EMG') %only .mat files with EMG (and skip flagEMG file)
            EMG{i}.EMGData=Datastr.EMG.EMGData;
            EMG{i}.EMGFrameRate=Datastr.EMG.EMGFrameRate;
            EMG{i}.EMGDataLabel=Datastr.EMG.EMGDataLabel;
            iDx=[iDx i];
            clear Datastr
        else
            warning(['No EMG data in file ' Datastr.Info.Trial '. Skipping.']);
        end
        end
    end
    for i=1:length(iDx)
        c=iDx(i);
        trialsList{i}=matList{c};
    end
    fprintf('\n#### Loading trial/s to get max EMG values... ####\n')
    [List,~]= uigetfile('*.mat','Select MVC trials/trials to get max EMG values','MultiSelect','on');
    
    if ~iscell(List) %Keep MaxEmgTrialsList as a cell, even if it is only one file
        MaxEmgTrialsList{1}=List;
    else
        MaxEmgTrialsList=List;
    end
    for iMax=1:size(MaxEmgTrialsList,2)
        load(MaxEmgTrialsList{iMax});
        EMGForMax{iMax}.EMGData=Datastr.EMG.EMGData;
        EMGForMax{iMax}.EMGDataLabel=Datastr.EMG.EMGDataLabel;
    end
    for k=1:length(trialsList)
        c=iDx(k);
        EMGselectionIndexes{k}=findIndexes(EMG{c}.EMGDataLabel,EMGs_Input);
        EMGsSelected{k}=EMG{c}.EMGData(:,EMGselectionIndexes{k});
    end
    
    
    for k=1:length(MaxEmgTrialsList)
        EMGselectionIndexesForMax{k}=findIndexes(EMGForMax{k}.EMGDataLabel,EMGs_Input);
        EMGsSelectedForMax{k}=EMGForMax{k}.EMGData(:,EMGselectionIndexesForMax{k});
    end
    
    
    
    AnalogFrameRate=EMG{iDx(1)}.EMGFrameRate; %Assuming all trials have same EMG rate
    %% ------------------------------------------------------------------------
    %                       EMG FILTERING: ENVELOPE
    %--------------------------------------------------------------------------
    fprintf('\n#### Filtering trials to be normalized... ####\n')
    %fcut for EMG assumed fixed (6Hz)
    EMGsEnvelope=EMGFiltering(EMGsSelected,AnalogFrameRate);
    
    EMGsEnvelopeForMax=EMGFiltering(EMGsSelectedForMax,AnalogFrameRate);
    
    %% ------------------------------------------------------------------------
    %                      EMG ANALYSIS WINDOW SELECTION
    %--------------------------------------------------------------------------
  
    %The analysis window for the trials is the whole trial now. Data is cropped in another module in BioMechPro:
    AnalysisWindow=[];
    EMGOffset=[];
    [EMGsFiltered,~]=selectionData_modified(EMGsEnvelope,AnalysisWindow,AnalogFrameRate,EMGOffset,0);

    
    %% ------------------------------------------------------------------------
    %                        COMPUTE MAX EMG VALUES
    %--------------------------------------------------------------------------
    fprintf('\n#### Computing max emg values... ####\n')
    [MaxEMG_aframes, numMaxEMG_trials,MaxEMGvalues]=computeMaxEMGvalues(EMGsEnvelopeForMax);
    
    sMaxEMG_trials=MaxEmgTrialsList(numMaxEMG_trials);
    MaxEMG_time=MaxEMG_aframes/AnalogFrameRate;
    
    %print maxemg.txt
    printMaxEMGvalues(maxemgPath, EMGs_Input, MaxEMGvalues, sMaxEMG_trials, MaxEMG_time);
    
    disp('Printed maxemg.txt')
    fprintf('\n#### Normalizing according to maxemg.txt... ####\n')
    %% ------------------------------------------------------------------------
    %                            NORMALIZE EMG
    %--------------------------------------------------------------------------
    if setOutliersTo1 % if EMG surpasses the maximum value found in MVC, it is set to 1
        NormEMG=normalizeEMG_mod(EMGsFiltered,MaxEMGvalues);
    else % values surpassing the max emg value will have values higher than 1
    NormEMG=normalizeEMG(EMGsFiltered,MaxEMGvalues);
    end
    %% ------------------------------------------------------------------------
    %                          SAVING and PLOTTING
    %--------------------------------------------------------------------------
    
    fprintf('\n#### Saving in Datastr... This might take a while ####\n')
    % Store data
    cd(Datastr.Info.SubjRoot)
    for i=1:length(trialsList)
        load(trialsList{i});
        Datastr.EMG.filtEMGData = EMGsFiltered{i};
        Datastr.EMG.normEMGData = NormEMG{i};
        %Datastr.EMG.filtEMGtime = EMGtime{i};
        save(trialsList{i},'Datastr')
    end
    disp('(flagEMG set to 1)')
end
cd(Datastr.Info.SubjRoot) %Make sure it is in subject folder
if ~flagEMG
    flagEMG=1;
    save('flagEMG.mat','flagEMG')
    % The first trial was processsed ((end-6:1:end-4)) assumes a filename with XXX.mat
    % with XXX being the trial number
    fprintf(['Trial ' matList{1}(end-6:1:end-4) ' processed (flagEMG=1).\n ']) 
end

Datastr=[]; %Datastr was previously saved for each trial. No need to do it here.
end