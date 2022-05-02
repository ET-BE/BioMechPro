function [Datastr] = S22_aveGait(Datastr, figurePath, EMGNorType, plotGen)
% gBMPDynUI figurePath=1; EMGNorType=1; plotGen=1;
%
% plot the joint motion data of two systems, after the sychorinization
% INPUT)
% - Datastr, the data structure with at least the fields:
%          .Resample.Sych.IKAngData
%          .Resample.Sych.IDTrqData
%          .Resample.Sych.IMUAngData
%          .Resample.Sych.IDAngData_Portable
%          .Resample.Sych.ForcePlateGRFData
%          .Resample.Sych.InsoleGRFData
%
% - trial_name, the trial names that will excute this average function

    %% check variables
    if ~ isfield(Datastr.Resample.Sych, 'IKAngData')
        warning('No inverse kinematics data for averaging')
    end

    if ~ isfield(Datastr.Resample.Sych, 'IMUAngData')
        warning('No IMU angle data for averaging')
    end

    if ~ isfield(Datastr.Resample.Sych, 'IDTrqData')
        warning('No inverse dynamics data for averaging')
    end

    if ~ isfield(Datastr.Resample.Sych, 'IDTrqData_Portable')
        warning('No inverse dynamics data of portable system for averaging')
    end

    if ~ isfield(Datastr.Resample.Sych, 'ForcePlateGRFData')
        warning('No force plate force data for averaging')
    end

    if ~ isfield(Datastr.Resample.Sych, 'InsoleGRFData')
        warning('No insole force data for averaging')
    end
    
    if ~ isfield(Datastr.Resample.Sych, 'ForcePlateGRFDataInCalcn')
        warning('No GRF data in Calcn for averaging')
    end
    
    %% average the gait cycle
    
    if strcmp(Datastr.Info.Trial(1:3), 'run') || strcmp(Datastr.Info.Trial(1:4), 'walk')
       % in walking and running trials, using Fy to detect gait cycles
       hsMatrix_l = Datastr.Resample.Sych.Average.hsMatrix_left;
       hsMatrix_r = Datastr.Resample.Sych.Average.hsMatrix_right;
    else
        hsMatrix_l = Datastr.Resample.Sych.Average.hsMatrix_right;
        hsMatrix_r = Datastr.Resample.Sych.Average.hsMatrix_right;
    end
       
%        speed = str2double(Datastr.Info.Trial(end-1:end))/36;

       % generate the figure save folder
       figSavePath = strcat(Datastr.Info.SubjRoot, '\', figurePath, Datastr.Info.Trial);
       if ~exist(figSavePath, 'dir')
           mkdir(figSavePath)
       end

       % calculate IKAngData average and std
       if isfield(Datastr.Resample.Sych, 'IKAngData')
           [interpVar_l, aveVar_l, stdVar_l, interpVar_r, aveVar_r, stdVar_r] = ...
                getGaitAve_noCheck(hsMatrix_l, hsMatrix_r, Datastr.Resample.Sych.IKAngData,...
                Datastr.Resample.FrameRate, 'IKAngData', figSavePath, plotGen);

            Datastr.Resample.Sych.Average.IKAngData.mat_l = interpVar_l;
            Datastr.Resample.Sych.Average.IKAngData.ave_l = aveVar_l;
            Datastr.Resample.Sych.Average.IKAngData.std_l = stdVar_l;
            Datastr.Resample.Sych.Average.IKAngData.mat_r = interpVar_r;
            Datastr.Resample.Sych.Average.IKAngData.ave_r = aveVar_r;
            Datastr.Resample.Sych.Average.IKAngData.std_r = stdVar_r;
       end

       if isfield(Datastr.Resample.Sych, 'IMUAngData')
           % calculate IMUAngData average and std
           [interpVar_l, aveVar_l, stdVar_l, interpVar_r, aveVar_r, stdVar_r] = ...
                getGaitAve_noCheck(hsMatrix_l, hsMatrix_r, Datastr.Resample.Sych.IMUAngData,...
                Datastr.Resample.FrameRate, 'IMUAngData', figSavePath, plotGen);

            Datastr.Resample.Sych.Average.IMUAngData.mat_l = interpVar_l;
            Datastr.Resample.Sych.Average.IMUAngData.ave_l = aveVar_l;
            Datastr.Resample.Sych.Average.IMUAngData.std_l = stdVar_l;
            Datastr.Resample.Sych.Average.IMUAngData.mat_r = interpVar_r;
            Datastr.Resample.Sych.Average.IMUAngData.ave_r = aveVar_r;
            Datastr.Resample.Sych.Average.IMUAngData.std_r = stdVar_r;
       end
       
       if isfield(Datastr.Resample.Sych, 'IDTrqData')
           % calculate IDTrqData average and std
           [interpVar_l, aveVar_l, stdVar_l, interpVar_r, aveVar_r, stdVar_r] = ...
                getGaitAve_noCheck(hsMatrix_l, hsMatrix_r, Datastr.Resample.Sych.IDTrqData,...
                Datastr.Resample.FrameRate, 'IDTrqData', figSavePath, plotGen);

            Datastr.Resample.Sych.Average.IDTrqData.mat_l = interpVar_l;
            Datastr.Resample.Sych.Average.IDTrqData.ave_l = aveVar_l;
            Datastr.Resample.Sych.Average.IDTrqData.std_l = stdVar_l;
            Datastr.Resample.Sych.Average.IDTrqData.mat_r = interpVar_r;
            Datastr.Resample.Sych.Average.IDTrqData.ave_r = aveVar_r;
            Datastr.Resample.Sych.Average.IDTrqData.std_r = stdVar_r;
       end
       
       if isfield(Datastr.Resample.Sych, 'IDTrqData_Portable')
           % calculate IDTrqData_Portable average and std
           [interpVar_l, aveVar_l, stdVar_l, interpVar_r, aveVar_r, stdVar_r] = ...
                getGaitAve_noCheck(hsMatrix_l, hsMatrix_r, Datastr.Resample.Sych.IDTrqData_Portable,...
                Datastr.Resample.FrameRate, 'IDTrqData_Portable', figSavePath, plotGen);

            Datastr.Resample.Sych.Average.IDTrqData_Portable.mat_l = interpVar_l;
            Datastr.Resample.Sych.Average.IDTrqData_Portable.ave_l = aveVar_l;
            Datastr.Resample.Sych.Average.IDTrqData_Portable.std_l = stdVar_l;
            Datastr.Resample.Sych.Average.IDTrqData_Portable.mat_r = interpVar_r;
            Datastr.Resample.Sych.Average.IDTrqData_Portable.ave_r = aveVar_r;
            Datastr.Resample.Sych.Average.IDTrqData_Portable.std_r = stdVar_r;
       end

       if isfield(Datastr.Resample.Sych, 'ForcePlateGRFData')
            % calculate ForcePlateGRFData average and std
           [interpVar_l, aveVar_l, stdVar_l, interpVar_r, aveVar_r, stdVar_r] = ...
                getGaitAve_noCheck(hsMatrix_l, hsMatrix_r, Datastr.Resample.Sych.ForcePlateGRFData,...
                Datastr.Resample.FrameRate, 'ForcePlateGRFData', figSavePath, plotGen);

            Datastr.Resample.Sych.Average.ForcePlateGRFData.mat_l = interpVar_l;
            Datastr.Resample.Sych.Average.ForcePlateGRFData.ave_l = aveVar_l;
            Datastr.Resample.Sych.Average.ForcePlateGRFData.std_l = stdVar_l;
            Datastr.Resample.Sych.Average.ForcePlateGRFData.mat_r = interpVar_r;
            Datastr.Resample.Sych.Average.ForcePlateGRFData.ave_r = aveVar_r;
            Datastr.Resample.Sych.Average.ForcePlateGRFData.std_r = stdVar_r;
       end

       if isfield(Datastr.Resample.Sych, 'InsoleGRFData')
            % calculate InsoleGRFData average and std
           [interpVar_l, aveVar_l, stdVar_l, interpVar_r, aveVar_r, stdVar_r] = ...
                getGaitAve_noCheck(hsMatrix_l, hsMatrix_r, Datastr.Resample.Sych.InsoleGRFData,...
                Datastr.Resample.FrameRate, 'InsoleGRFData', figSavePath, plotGen);

            Datastr.Resample.Sych.Average.InsoleGRFData.mat_l = interpVar_l;
            Datastr.Resample.Sych.Average.InsoleGRFData.ave_l = aveVar_l;
            Datastr.Resample.Sych.Average.InsoleGRFData.std_l = stdVar_l;
            Datastr.Resample.Sych.Average.InsoleGRFData.mat_r = interpVar_r;
            Datastr.Resample.Sych.Average.InsoleGRFData.ave_r = aveVar_r;
            Datastr.Resample.Sych.Average.InsoleGRFData.std_r = stdVar_r;
       end
       
       if isfield(Datastr.Resample.Sych, 'ForcePlateGRFDataInCalcn')
           % calculate ForcePlateGRFData average and std
           [interpVar_l, aveVar_l, stdVar_l, interpVar_r, aveVar_r, stdVar_r] = ...
                getGaitAve_noCheck(hsMatrix_l, hsMatrix_r, Datastr.Resample.Sych.ForcePlateGRFDataInCalcn,...
                Datastr.Resample.FrameRate, 'ForcePlateGRFDataInCalcn', figSavePath, plotGen);

            Datastr.Resample.Sych.Average.ForcePlateGRFDataInCalcn.mat_l = interpVar_l;
            Datastr.Resample.Sych.Average.ForcePlateGRFDataInCalcn.ave_l = aveVar_l;
            Datastr.Resample.Sych.Average.ForcePlateGRFDataInCalcn.std_l = stdVar_l;
            Datastr.Resample.Sych.Average.ForcePlateGRFDataInCalcn.mat_r = interpVar_r;
            Datastr.Resample.Sych.Average.ForcePlateGRFDataInCalcn.ave_r = aveVar_r;
            Datastr.Resample.Sych.Average.ForcePlateGRFDataInCalcn.std_r = stdVar_r;
       end
       % calculate InsoleGRFData average and std
       if strcmpi(EMGNorType, 'dynMVC')
           dynMVCvalue = importdata([Datastr.Info.SubjRoot, '\dynMVCvalue.mat']);

           [interpVar_l, aveVar_l, stdVar_l, interpVar_r, aveVar_r, stdVar_r] = ...
            getGaitAve_noCheck(hsMatrix_l, hsMatrix_r, Datastr.Resample.Sych.EMG./dynMVCvalue,...
            Datastr.Resample.FrameRate, 'EMG', figSavePath, plotGen);
            Datastr.Resample.Sych.Average.EMGAvedynNorFlag = 1;
            Datastr.Resample.Sych.Average.EMGAveNorFlag = 0;
        
       elseif strcmpi(EMGNorType, 'MVC')
           MVCvalue = importdata([Datastr.Info.SubjRoot, '\MVCvalue.mat']);

           [interpVar_l, aveVar_l, stdVar_l, interpVar_r, aveVar_r, stdVar_r] = ...
            getGaitAve_noCheck(hsMatrix_l, hsMatrix_r, Datastr.Resample.Sych.EMG./MVCvalue,...
            Datastr.Resample.FrameRate, 'EMG', figSavePath, plotGen);
            Datastr.Resample.Sych.Average.EMGAvedynNorFlag = 0;
            Datastr.Resample.Sych.Average.EMGAveNorFlag = 1;
        
       else
           [interpVar_l, aveVar_l, stdVar_l, interpVar_r, aveVar_r, stdVar_r] = ...
            getGaitAve_noCheck(hsMatrix_l, hsMatrix_r, Datastr.Resample.Sych.EMG,...
            Datastr.Resample.FrameRate, 'EMG', figSavePath, plotGen);
            Datastr.Resample.Sych.Average.EMGAvedynNorFlag = 0;
            Datastr.Resample.Sych.Average.EMGAveNorFlag = 0;
       end

        Datastr.Resample.Sych.Average.EMG.mat_l = interpVar_l;
        Datastr.Resample.Sych.Average.EMG.ave_l = aveVar_l;
        Datastr.Resample.Sych.Average.EMG.std_l = stdVar_l;
        Datastr.Resample.Sych.Average.EMG.mat_r = interpVar_r;
        Datastr.Resample.Sych.Average.EMG.ave_r = aveVar_r;
        Datastr.Resample.Sych.Average.EMG.std_r = stdVar_r;
    
    %% plot the envelop of mean and std
    if strcmpi(plotGen, 'True') 
        fig_save_path = 'Figures\aveMeanStd';
        plotAveMeanStd(Datastr, fig_save_path);
    else
        fprintf('plot not generated, please put "True" in the "PlotGen" if you would like to have plot generated!/n')
    end

end



