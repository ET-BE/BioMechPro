%% Process 3
% Cut data into pieces and sort it
%
% NOTE : THIS IS AN EXPERIMENT SPECIFIC FILE

clear variables; close all; clc;

subjects = 1:10;  % corresponding with how many you have in getInfo

%% Cut data on perturbation events

for isubj = subjects

    % Clean up some variables
    clear pieceMarker pieceEvent pieceForce pieceEMG pieceOther trialSet
    
    % Get required information
    getInfo;
    
    % Get trial info
    switch isubj
        case 1
            trialSet(1,:) = [101:103]; tableSet(1) = 1; % Slow ML
            trialSet(2,:) = [104:106]; tableSet(2) = 1; % Fast ML
            trialSet(3,:) = [107:109]; tableSet(3) = 2; % Slow AP
            trialSet(4,:) = [110:112]; tableSet(4) = 2; % Fast AP
            subjTblFile = '141031PP1randTbl';
        case 2
            trialSet(1,1:5) = [101:105]; tableSet(1) = 1;
            trialSet(2,1:3) = [106:108]; tableSet(2) = 1;
            trialSet(3,1:3) = [109:111]; tableSet(3) = 2; 
            trialSet(4,1:3) = [112:114]; tableSet(4) = 2; 
            subjTblFile = '141106PP2randTbl';
        case 3
            trialSet(1,:) = [101:103]; tableSet(1) = 1; 
            trialSet(2,:) = [104:106]; tableSet(2) = 1;
            trialSet(3,:) = [107:109]; tableSet(3) = 2;
            trialSet(4,:) = [110:112]; tableSet(4) = 2;
            subjTblFile = '141106PP3randTbl';
        case 4
            trialSet(1,:) = [101:103]; tableSet(1) = 1; 
            trialSet(2,:) = [104:106]; tableSet(2) = 1;
            trialSet(3,:) = [107:109]; tableSet(3) = 2; 
            trialSet(4,:) = [110:112]; tableSet(4) = 2;
            subjTblFile = '141107PP4randTbl';
        case 5
            trialSet(1,:) = [101:103]; tableSet(1) = 1;
            trialSet(2,:) = [104:106]; tableSet(2) = 1;
            trialSet(3,:) = [107:109]; tableSet(3) = 2;
            trialSet(4,:) = [110:112]; tableSet(4) = 2;
            subjTblFile = '141107PP5randTbl';
        case 6
            trialSet(1,:) = [101:103]; tableSet(1) = 1;
            trialSet(2,:) = [104:106]; tableSet(2) = 1;
            trialSet(3,:) = [107:109]; tableSet(3) = 2;
            trialSet(4,:) = [110:112]; tableSet(4) = 2;
            subjTblFile = '141119PP6randTbl';
        case 7
            trialSet(1,:) = [101:103]; tableSet(1) = 1;
            trialSet(2,:) = [104:106]; tableSet(2) = 1; 
            trialSet(3,:) = [107:109]; tableSet(3) = 2;
            trialSet(4,:) = [110:112]; tableSet(4) = 2;
            subjTblFile = '141125PP7randTbl';
        case 8
            trialSet(1,:) = [101:103]; tableSet(1) = 1;
            trialSet(2,:) = [104:106]; tableSet(2) = 1;
            trialSet(3,:) = [107:109]; tableSet(3) = 2;
            trialSet(4,:) = [110:112]; tableSet(4) = 2;
            subjTblFile = '141125PP8randTbl';
        case 9
            trialSet(1,1:3) = [101:103]; tableSet(1) = 1;
            trialSet(2,1:3) = [104:106]; tableSet(2) = 1;
            trialSet(3,1:3) = [107:109]; tableSet(3) = 2;
            trialSet(4,1:4) = [110:113]; tableSet(4) = 2;
            subjTblFile = '141126PP9randTbl';
        case 10
            trialSet(1,:) = [101:103]; tableSet(1) = 1;
            trialSet(2,:) = [104:106]; tableSet(2) = 1;
            trialSet(3,:) = [107:109]; tableSet(3) = 2;
            trialSet(4,:) = [110:112]; tableSet(4) = 2;
            subjTblFile = '141126PP10randTbl';            
        otherwise
            return
    end
    
    
    for iSet = 1:size(trialSet,1)
%     for iSet = [1 3]; % Slow only!
        % Load perturbation table and classification table 
        % classTbl is in Newton, convert to bodyprct , columns are 1:mag1, 2:mag2 3:del 4:dur 5:typenum
        load([subjroot '\XPC\' subjTblFile num2str(tableSet(iSet)) '.mat']);
        classTbl(:,[1,2]) = classTbl(:,[1 2]) ./ (subjmass*9.81);
        
        % Instantiate data containers (one for each field)
        if ~exist('pieceMarker','var')
            for iSet2 = 1:size(trialSet,1)
                pieceMarker.Set(iSet2).Base.Data = struct();
                pieceEvent.Set(iSet2).Base.Data = struct();
                pieceForce.Set(iSet2).Base.Data = struct();
                pieceEMG.Set(iSet2).Base.Data = struct();
                pieceOther.Set(iSet2).Base.Data = struct();
                
                pieceMarker.Set(iSet2).Pert(size(classTbl,1)).Data = struct();
                pieceEvent.Set(iSet2).Pert(size(classTbl,1)).Data = struct();
                pieceForce.Set(iSet2).Pert(size(classTbl,1)).Data = struct();
                pieceEMG.Set(iSet2).Pert(size(classTbl,1)).Data = struct();
                pieceOther.Set(iSet2).Pert(size(classTbl,1)).Data = struct();
            end
        end
        
        % Loop through trials for processing
        n1=0; n2=0; n3=0; n4=0; n5=0; n6=0; n7=0; n8=0;
        for itrial = setdiff(trialSet(iSet,:),0) % remove any zeros (can occur due filling: unequal trials/set)

            % Load trial
            backslsh = strfind(subjroot,'\');
            load([subjroot '\' subjroot(backslsh(end)+1:end) '_' num2str(itrial) '.mat'])

            % Get sample frequencies and offsets
            fs_mark = Datastr.Marker.MarkerFrameRate;
            fs_other = Datastr.Other.OtherFrameRate;
            fs_rate = fs_other / fs_mark;
            markerSyncIdx = Datastr.Marker.MarkerSyncIdx;
            otherSyncIdx = Datastr.Other.OtherSyncIdx; % = force sync = emg sync

            % Get events from perturbations
            % Note: + 2, because 1 sample delay with input trq and 1 sample correction for diff
            odata = Datastr.Other.OtherData;
            pertOnset = find( diff( odata(:,2) ) ) + 2;

            % Find indices of perturbation end / find perturbation magnitude (BW ratio)
            % NOTE : positive is a pull to the right
            if ~isempty(pertOnset)

                pertEnd = zeros(length(pertOnset),1);
                pertMag1 = zeros(length(pertOnset),1); % Magnitude vector
                pertMag2 = zeros(length(pertOnset),1); % Magnitude vector

                for i = 1:length(pertOnset) % For all perturbations
                    pertMag1(i,1) = odata(pertOnset(i),3) ./ (0.3 * subjmass * 9.81); % Input trq 1 to bw fraction
                    pertMag2(i,1) = odata(pertOnset(i),4) ./ (0.3 * subjmass * 9.81); % Input trq 2 to bw fraction

                    if pertMag1(i,1) ~= 0       % Perturbation on motor 1

                        if i ~= length(pertOnset)
                            pertEnd(i,1) = pertOnset(i) -1 + find( odata(pertOnset(i):pertOnset(i+1)-1 , 3) == 0 , 1 , 'first'); % Find where the input torque drops to 0
                        else % if there is no subsequent pert
                            pertEnd(i,1) = pertOnset(i) -1 + find( odata(pertOnset(i):end , 3) == 0 , 1 , 'first');
                        end

                    elseif pertMag2(i,1) ~= 0   % Perturbation on motor 2

                        if i ~= length(pertOnset)
                            pertEnd(i,1) = pertOnset(i) -1 + find( odata(pertOnset(i):pertOnset(i+1)-1 , 4) == 0 , 1 , 'first'); % Find where the input torque drops to 0
                        else % if there is no subsequent pert
                            pertEnd(i,1) = pertOnset(i) -1 + find( odata(pertOnset(i):end , 4) == 0 , 1 , 'first');
                        end

                    end

                end
                pertDur = (pertEnd-pertOnset)./fs_other;

                % Convert pertOnset to onset index suitable for marker data
%                     pertOnsetM = round( (pertOnset - otherSyncIdx + 1) ./ fs_rate )  + markerSyncIdx - 1;

                % Generate start and stop indices for cutting
                % OPTION 1: 3 second range with base = pertOnset
                baseIdx = pertOnset;
                idxStartM = round( (baseIdx - otherSyncIdx(1) + 1) ./ fs_rate )  + markerSyncIdx(1) - 1 - 0.5.*fs_mark + 1;
                idxStopM = round( (baseIdx - otherSyncIdx(1) + 1) ./ fs_rate )  + markerSyncIdx(1) - 1 + 2.5*fs_mark;
                idxStartO = baseIdx - 0.5.*fs_other + 1;
                idxStopO = baseIdx + 2.5*fs_other;

                % OPTION 2: two single samples, one at a given time after pertOnset, one at SSLend after pertOnset (first RHS)
%                     baseIdx = pertOnset + fs_other./10; % 100 ms after pert onset
%                     idxStartM = round( (baseIdx - otherSyncIdx + 1) ./ fs_rate )  + markerSyncIdx - 1;
%                     idxStopM = round( (baseIdx - otherSyncIdx + 1) ./ fs_rate )  + markerSyncIdx - 1;
%                     idxStartO = baseIdx;
%                     idxStopO = baseIdx;


                % Classify perturbations on their magnitude and duration, then CUT THE DATA
                for i = 1:length(pertDur)

                    % Classify perturbation on duration and magnitude
                    [~,ipert] = min( abs(classTbl(:,1) - pertMag1(i)) + abs(classTbl(:,2) - pertMag2(i)) + abs(classTbl(:,4) - pertDur(i)) );

                    % Count perturbations
                    switch ipert
                        case 1
                            n1 = n1+1;
                        case 2
                            n2 = n2+1;
                        case 3
                            n3 = n3+1;
                        case 4
                            n4 = n4+1;
                        case 5
                            n5 = n5+1;
                        case 6
                            n6 = n6+1;
                        case 7 
                            n7 = n7+1;
                        case 8
                            n8 = n8+1;
                    end
                    
                    % Cut and store the data you want
                    dataRangeM = idxStartM(i):idxStopM(i);
                    dataRangeO = idxStartO(i):idxStopO(i);

                    if ~isfield(pieceMarker.Set(iSet).Pert(ipert).Data,'MarkerFrameRate')
                        % Marker
                        pieceMarker.Set(iSet).Pert(ipert).Data.MarkerFrameRate = Datastr.Marker.MarkerFrameRate;
                        pieceMarker.Set(iSet).Pert(ipert).Data.COMData = Datastr.Marker.COMData(dataRangeM,:,:);
                        pieceMarker.Set(iSet).Pert(ipert).Data.COMDataD = Datastr.Marker.COMDataD(dataRangeM,:,:);
                        pieceMarker.Set(iSet).Pert(ipert).Data.COMDataLabel = Datastr.Marker.COMDataLabel;
%                         pieceMarker.Set(iSet).Pert(ipert).Data.ClusterAngleE = Datastr.Marker.ClusterAngleE(dataRangeM,:,:);
%                         pieceMarker.Set(iSet).Pert(ipert).Data.ClusterAngleED = Datastr.Marker.ClusterAngleED(dataRangeM,:,:);
%                         pieceMarker.Set(iSet).Pert(ipert).Data.ClusterAngleEDD = Datastr.Marker.ClusterAngleEDD(dataRangeM,:,:);
                        try
                            pieceMarker.Set(iSet).Pert(ipert).Data.JointAngData = Datastr.Marker.JointAngData(dataRangeM,:,:);
                            pieceMarker.Set(iSet).Pert(ipert).Data.JointAngDataD = Datastr.Marker.JointAngDataD(dataRangeM,:,:);
                            pieceMarker.Set(iSet).Pert(ipert).Data.JointAngDataLabel = Datastr.Marker.JointAngDataLabel;
                            pieceMarker.Set(iSet).Pert(ipert).Data.JointTrqData = Datastr.Marker.JointTrqData(dataRangeM,:,:);
                            pieceMarker.Set(iSet).Pert(ipert).Data.JointTrqDataLabel = Datastr.Marker.JointTrqDataLabel;
                        catch
                            pieceMarker.Set(iSet).Pert(ipert).Data.JointAngData = zeros([300 9 3 2]);
                            pieceMarker.Set(iSet).Pert(ipert).Data.JointAngDataD = zeros([300 9 3 2]);
                            pieceMarker.Set(iSet).Pert(ipert).Data.JointTrqData = zeros([300 9 3 2]);
                        end
                        
                        % Event
                        pieceEvent.Set(iSet).Pert(ipert).Data.GaitPhaseM = Datastr.Event.GaitPhaseM(dataRangeM,:);
                        pieceEvent.Set(iSet).Pert(ipert).Data.GaitPhaseMLabel = Datastr.Event.GaitPhaseMLabel;
                        pieceEvent.Set(iSet).Pert(ipert).Data.GaitPhaseF = Datastr.Event.GaitPhaseF(dataRangeO,:);
                        pieceEvent.Set(iSet).Pert(ipert).Data.GaitPhaseFLabel = Datastr.Event.GaitPhaseFLabel;

                        % Force
                        pieceForce.Set(iSet).Pert(ipert).Data.ForceData = Datastr.Force.ForceData(dataRangeO,:);
                        pieceForce.Set(iSet).Pert(ipert).Data.ForceDataLabel = Datastr.Force.ForceDataLabel;
                        pieceForce.Set(iSet).Pert(ipert).Data.ForceFrameRate = Datastr.Force.ForceFrameRate;

                        % EMG
                        pieceEMG.Set(iSet).Pert(ipert).Data.EMGData = Datastr.EMG.EMGData(dataRangeO,:);
                        pieceEMG.Set(iSet).Pert(ipert).Data.EMGDataLabel = Datastr.EMG.EMGDataLabel;
                        pieceEMG.Set(iSet).Pert(ipert).Data.EMGFrameRate = Datastr.EMG.EMGFrameRate;

                        % Other
                        pieceOther.Set(iSet).Pert(ipert).Data.OtherData = Datastr.Other.OtherData(dataRangeO,:);
                        pieceOther.Set(iSet).Pert(ipert).Data.OtherDataLabel = Datastr.Other.OtherDataLabel;
                        pieceOther.Set(iSet).Pert(ipert).Data.OtherFrameRate = Datastr.Other.OtherFrameRate;

                    else
                        % Marker
                        pieceMarker.Set(iSet).Pert(ipert).Data.COMData(:,:,:,end+1) = Datastr.Marker.COMData(dataRangeM,:,:);
                        pieceMarker.Set(iSet).Pert(ipert).Data.COMDataD(:,:,:,end+1) = Datastr.Marker.COMDataD(dataRangeM,:,:);
%                         pieceMarker.Set(iSet).Pert(ipert).Data.ClusterAngleE(:,:,:,end+1) = Datastr.Marker.ClusterAngleE(dataRangeM,:,:);
%                         pieceMarker.Set(iSet).Pert(ipert).Data.ClusterAngleED(:,:,:,end+1) = Datastr.Marker.ClusterAngleED(dataRangeM,:,:);
%                         pieceMarker.Set(iSet).Pert(ipert).Data.ClusterAngleEDD(:,:,:,end+1) = Datastr.Marker.ClusterAngleEDD(dataRangeM,:,:);
                        try
                            pieceMarker.Set(iSet).Pert(ipert).Data.JointAngData(:,:,:,end+1) = Datastr.Marker.JointAngData(dataRangeM,:,:);
                            pieceMarker.Set(iSet).Pert(ipert).Data.JointAngDataD(:,:,:,end+1) = Datastr.Marker.JointAngDataD(dataRangeM,:,:);
                            pieceMarker.Set(iSet).Pert(ipert).Data.JointTrqData(:,:,:,end+1) = Datastr.Marker.JointTrqData(dataRangeM,:,:);
                        catch
%                             disp(['Empty pieceMarker ' num2str(isubj)]);
                        end
                        
                        % Event
                        pieceEvent.Set(iSet).Pert(ipert).Data.GaitPhaseM(:,:,end+1) = Datastr.Event.GaitPhaseM(dataRangeM,:);
                        pieceEvent.Set(iSet).Pert(ipert).Data.GaitPhaseF(:,:,end+1) = Datastr.Event.GaitPhaseF(dataRangeO,:);
                        
                        % Force
                        pieceForce.Set(iSet).Pert(ipert).Data.ForceData(:,:,end+1) = Datastr.Force.ForceData(dataRangeO,:);

                        % EMG
                        pieceEMG.Set(iSet).Pert(ipert).Data.EMGData(:,:,end+1) = Datastr.EMG.EMGData(dataRangeO,:);

                        % Other
                        pieceOther.Set(iSet).Pert(ipert).Data.OtherData(:,:,end+1) = Datastr.Other.OtherData(dataRangeO,:);
                    end
                end

            else % IS BASELINE TRIAL

                % Generate start and stop indices for cutting
                % OPTION 1: 3 second range with base = SSLstart
%                 baseIdx = Datastr.Event.GaitPhaseMIdx.SSL(:,1); % SSLstart
%                 idxStartM = baseIdx - 0.5.*fs_mark + 1; 
%                 idxStopM = baseIdx + 2.5.*fs_mark;
%                 idxStartO = otherSyncIdx - 1 + (baseIdx - markerSyncIdx).*10 - 0.5.*fs_other + 1;
%                 idxStopO = otherSyncIdx - 1 + (baseIdx - markerSyncIdx).*10 + 2.5.*fs_other;

                % OPTION 2: 3 second range with base = SSLstart from Force data
                baseIdx = Datastr.Event.GaitPhaseFIdx.SSL(:,1); % SSLstart
                corrIdx = 0; % 20
                idxStartM = round( (baseIdx - otherSyncIdx(1) + 1 - corrIdx) ./ fs_rate )  + markerSyncIdx(1) - 1 - 0.5.*fs_mark + 1;
                idxStopM = round( (baseIdx - otherSyncIdx(1) + 1 - corrIdx) ./ fs_rate )  + markerSyncIdx(1) - 1 + 2.5*fs_mark;
                idxStartO = baseIdx - 0.5.*fs_other + 1 - corrIdx;
                idxStopO = baseIdx + 2.5.*fs_other - corrIdx;
                
                % Don't take the first and last 10 steps
                iskip = 10;

                % Cut and store the data you want
                for i = iskip:length(idxStartM)-iskip

                    dataRangeM = idxStartM(i):idxStopM(i);
                    dataRangeO = idxStartO(i):idxStopO(i);

                    if ~isfield(pieceMarker.Set(iSet).Base.Data,'MarkerFrameRate')
                        % Marker
                        pieceMarker.Set(iSet).Base.Data.MarkerFrameRate = Datastr.Marker.MarkerFrameRate;
                        pieceMarker.Set(iSet).Base.Data.COMData = Datastr.Marker.COMData(dataRangeM,:,:);
                        pieceMarker.Set(iSet).Base.Data.COMDataD = Datastr.Marker.COMDataD(dataRangeM,:,:);
                        pieceMarker.Set(iSet).Base.Data.COMDataLabel = Datastr.Marker.COMDataLabel;
%                         pieceMarker.Set(iSet).Base.Data.ClusterAngleE = Datastr.Marker.ClusterAngleE(dataRangeM,:,:);
%                         pieceMarker.Set(iSet).Base.Data.ClusterAngleED = Datastr.Marker.ClusterAngleED(dataRangeM,:,:);
%                         pieceMarker.Set(iSet).Base.Data.ClusterAngleEDD = Datastr.Marker.ClusterAngleEDD(dataRangeM,:,:);
                        try
                            pieceMarker.Set(iSet).Base.Data.JointAngData = Datastr.Marker.JointAngData(dataRangeM,:,:);
                            pieceMarker.Set(iSet).Base.Data.JointAngDataD = Datastr.Marker.JointAngDataD(dataRangeM,:,:);
                            pieceMarker.Set(iSet).Base.Data.JointAngDataLabel = Datastr.Marker.JointAngDataLabel;
                            pieceMarker.Set(iSet).Base.Data.JointTrqData = Datastr.Marker.JointTrqData(dataRangeM,:,:);
                            pieceMarker.Set(iSet).Base.Data.JointTrqDataLabel = Datastr.Marker.JointTrqDataLabel;
                        catch
                            pieceMarker.Set(iSet).Base.Data.JointAngData = zeros([300 9 3 2]);
                            pieceMarker.Set(iSet).Base.Data.JointAngDataD = zeros([300 9 3 2]);
                            pieceMarker.Set(iSet).Base.Data.JointTrqData = zeros([300 9 3 2]);
                        end
                        
                        % Event
                        pieceEvent.Set(iSet).Base.Data.GaitPhaseM = Datastr.Event.GaitPhaseM(dataRangeM,:);
                        pieceEvent.Set(iSet).Base.Data.GaitPhaseMLabel = Datastr.Event.GaitPhaseMLabel;
                        pieceEvent.Set(iSet).Base.Data.GaitPhaseF = Datastr.Event.GaitPhaseF(dataRangeO,:);
                        pieceEvent.Set(iSet).Base.Data.GaitPhaseFLabel = Datastr.Event.GaitPhaseFLabel;
                        
                        % Force
                        pieceForce.Set(iSet).Base.Data.ForceData = Datastr.Force.ForceData(dataRangeO,:);
                        pieceForce.Set(iSet).Base.Data.ForceDataLabel = Datastr.Force.ForceDataLabel;
                        pieceForce.Set(iSet).Base.Data.ForceFrameRate = Datastr.Force.ForceFrameRate;
                        pieceForce.Set(iSet).Base.Data.ForceSyncIdx = Datastr.Force.ForceSyncIdx;
                        
                        % EMG
                        pieceEMG.Set(iSet).Base.Data.EMGData = Datastr.EMG.EMGData(dataRangeO,:);
                        pieceEMG.Set(iSet).Base.Data.EMGDataLabel = Datastr.EMG.EMGDataLabel;
                        pieceEMG.Set(iSet).Base.Data.EMGFrameRate = Datastr.EMG.EMGFrameRate;

                        % Other
                        pieceOther.Set(iSet).Base.Data.OtherData = Datastr.Other.OtherData(dataRangeO,:);
                        pieceOther.Set(iSet).Base.Data.OtherDataLabel = Datastr.Other.OtherDataLabel;
                        pieceOther.Set(iSet).Base.Data.OtherFrameRate = Datastr.Other.OtherFrameRate;

                    else
                        % Marker
                        pieceMarker.Set(iSet).Base.Data.COMDataD(:,:,:,end+1) = Datastr.Marker.COMDataD(dataRangeM,:,:);
                        pieceMarker.Set(iSet).Base.Data.COMData(:,:,:,end+1) = Datastr.Marker.COMData(dataRangeM,:,:);
%                         pieceMarker.Set(iSet).Base.Data.ClusterAngleE(:,:,:,end+1) = Datastr.Marker.ClusterAngleE(dataRangeM,:,:);
%                         pieceMarker.Set(iSet).Base.Data.ClusterAngleED(:,:,:,end+1) = Datastr.Marker.ClusterAngleED(dataRangeM,:,:);
%                         pieceMarker.Set(iSet).Base.Data.ClusterAngleEDD(:,:,:,end+1) = Datastr.Marker.ClusterAngleEDD(dataRangeM,:,:);
                        try
                            pieceMarker.Set(iSet).Base.Data.JointAngData(:,:,:,end+1) = Datastr.Marker.JointAngData(dataRangeM,:,:);
                            pieceMarker.Set(iSet).Base.Data.JointAngDataD(:,:,:,end+1) = Datastr.Marker.JointAngDataD(dataRangeM,:,:);
                            pieceMarker.Set(iSet).Base.Data.JointTrqData(:,:,:,end+1) = Datastr.Marker.JointTrqData(dataRangeM,:,:);
                        catch
%                             disp(['Empty pieceMarker ' num2str(isubj)]);
                        end
                        
                        % Event
                        pieceEvent.Set(iSet).Base.Data.GaitPhaseM(:,:,end+1) = Datastr.Event.GaitPhaseM(dataRangeM,:);
                        pieceEvent.Set(iSet).Base.Data.GaitPhaseF(:,:,end+1) = Datastr.Event.GaitPhaseF(dataRangeO,:);

                        % Force
                        pieceForce.Set(iSet).Base.Data.ForceData(:,:,end+1) = Datastr.Force.ForceData(dataRangeO,:);

                        % EMG
                        pieceEMG.Set(iSet).Base.Data.EMGData(:,:,end+1) = Datastr.EMG.EMGData(dataRangeO,:);

                        % Other
                        pieceOther.Set(iSet).Base.Data.OtherData(:,:,end+1) = Datastr.Other.OtherData(dataRangeO,:);
                    end

                end

            end

            disp(['Trial ' num2str(itrial) ' done']);
        end
            
    end

    disp([n1 n2 n3 n4 n5 n6 n7 n8]);
    
    % Save cut data into A NEW FILE 
    disp('Saving...might take a while....');
    
    save([subjroot '\pieceMarker'] , 'pieceMarker');  % The second time the subjfolder is used as filename
    save([subjroot '\pieceEvent'] , 'pieceEvent');  % The second time the subjfolder is used as filename
    save([subjroot '\pieceForce'] , 'pieceForce');  % The second time the subjfolder is used as filename
    save([subjroot '\pieceEMG'] , 'pieceEMG');  % The second time the subjfolder is used as filename
    save([subjroot '\pieceOther'] , 'pieceOther');  % The second time the subjfolder is used as filename
    
    disp(['isubj = ' num2str(isubj) ' done']);
    
end
