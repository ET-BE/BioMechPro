function plotAveMeanStd(Datastr, fig_save_path)
% gBMPDynUI trial_name=1; fig_save_path=1;
%
% plot the joint motion data of two systems, after the sychorinization
% INPUT)
% - Datastr, the data structure with at least the fields:
%          .Resample.IMU
%          .Resample.IKID
%          .Resample.Force
%          .Resample.Insole
%
% - MotionSychSign, the sign of doing motion sychronization or not
% - ForceSychSign, the sign of doing force sychronization or not

% Note: data will be plotted at multiple plots:
%   Fig11: joint motions 
%   Fig12: joint moments 
%   Fig13: vertical forces/CoP
%   Fig14: EMGs

%% check input

if isempty(fig_save_path)
    fig_save_path = '\Figures\aveMeanStd\';
end

% if strcmpi(trial_name, 'All')
%     plot_sign = 1;
% elseif any(strcmpi(Datastr.Info.Trial, trial_name))
%     plot_sign = 1;
% else
%     plot_sign = 0;
% end

    
if ~isfolder(strcat(Datastr.Info.SubjRoot, '\', fig_save_path))
    mkdir(strcat(Datastr.Info.SubjRoot, '\', fig_save_path))
end

joint_labels = ["LHip", "RHip", "LKnee", "RKnee", "LAnkle", "RAnkle"];
joint_angle_id = [15, 8, 18, 11, 19, 12];
joint_torque_id = [11, 8, 18, 17, 20, 19];

phase = 1:100;

figure()   % plot joint angles
for j = 1:6
    subplot(3, 2, j)
    if rem(j, 2) == 1
        errorfill(phase,...
             Datastr.Resample.Sych.Average.IKAngData.ave_l(:, joint_angle_id(j))',...
             Datastr.Resample.Sych.Average.IKAngData.std_l(:, joint_angle_id(j))', 'r-', 1.0);
        hold on;
        errorfill(phase,...
             Datastr.Resample.Sych.Average.IMUAngData.ave_l(:, joint_angle_id(j))',...
             Datastr.Resample.Sych.Average.IMUAngData.std_l(:, joint_angle_id(j))', 'b-', 0.6);
        hold off
    else
        errorfill(phase,...
             Datastr.Resample.Sych.Average.IKAngData.ave_r(:, joint_angle_id(j))',...
             Datastr.Resample.Sych.Average.IKAngData.std_r(:, joint_angle_id(j))', 'r-', 1.0);
        hold on;
        errorfill(phase,...
             Datastr.Resample.Sych.Average.IMUAngData.ave_r(:, joint_angle_id(j))',...
             Datastr.Resample.Sych.Average.IMUAngData.std_r(:, joint_angle_id(j))', 'b-', 0.6);
        hold off
    end
    title(joint_labels(j))
    if j == 5 || j == 6
        xlabel('Phase (%)')
    end
    if j == 1|| j == 3|| j == 5
        ylabel('Angle (deg.)')
    end
    if j == 2
        legend('Optical-std', 'Optical-mean','IMU-std', 'IMU-mean')
    end
    xlim([0, 101])
end
sgtitle(strcat('Joint Angle - ', Datastr.Info.Trial))
savefig(strcat(Datastr.Info.SubjRoot, '\', fig_save_path, '\JointAngle-', Datastr.Info.Trial, '.fig'))

figure()  % plot joint torques
for j = 1:6
    subplot(3, 2, j)
    if rem(j, 2) == 1
        errorfill(phase,...
             Datastr.Resample.Sych.Average.IDTrqData.ave_l(:, joint_torque_id(j))',...
             Datastr.Resample.Sych.Average.IDTrqData.std_l(:, joint_torque_id(j))', 'r-', 1.0);
        hold on;
        errorfill(phase,...
             Datastr.Resample.Sych.Average.IDTrqData_Portable.ave_l(:, joint_torque_id(j))',...
             Datastr.Resample.Sych.Average.IDTrqData_Portable.std_l(:, joint_torque_id(j))', 'b-', 0.6);
        hold off
    else
        errorfill(phase,...
             Datastr.Resample.Sych.Average.IDTrqData.ave_r(:, joint_torque_id(j))',...
             Datastr.Resample.Sych.Average.IDTrqData.std_r(:, joint_torque_id(j))', 'r-', 1.0);
        hold on;
        errorfill(phase,...
             Datastr.Resample.Sych.Average.IDTrqData_Portable.ave_r(:, joint_torque_id(j))',...
             Datastr.Resample.Sych.Average.IDTrqData_Portable.std_r(:, joint_torque_id(j))', 'b-', 0.6);
        hold off
    end

    title(joint_labels(j))
    if j == 5 || j == 6
        xlabel('Phase (%)')
    end
    if j == 1|| j == 3|| j == 5
        ylabel('Torque (Nm)')
    end
    if j == 2
        legend('Optical-std', 'Optical-mean', 'Portable-std', 'Portable-mean')
    end
    xlim([0, 101])
end
sgtitle(strcat('Joint Torque - ', Datastr.Info.Trial))
savefig(strcat(Datastr.Info.SubjRoot, '\', fig_save_path, '\JointTorque-', Datastr.Info.Trial, '.fig'))

figure()  % plot Fy and CoPs
subplot(2, 3, 1)
errorfill(phase,...
     Datastr.Resample.Sych.Average.ForcePlateGRFDataInCalcn.ave_l(:, 3)',...
     Datastr.Resample.Sych.Average.ForcePlateGRFDataInCalcn.std_l(:, 3)', 'r-', 1.0);
hold on;
errorfill(phase,...
     Datastr.Resample.Sych.Average.InsoleGRFData.ave_l(:, 3)',...
     Datastr.Resample.Sych.Average.InsoleGRFData.std_l(:, 3)', 'b-', 0.6);
hold off
xlim([0, 101])
title('L Fy')
ylabel('N')

subplot(2, 3, 2)
errorfill(phase,...
     Datastr.Resample.Sych.Average.ForcePlateGRFDataInCalcn.ave_l(:, 5)',...
     Datastr.Resample.Sych.Average.ForcePlateGRFDataInCalcn.std_l(:, 5)', 'r-', 1.0);
hold on
errorfill(phase,...
     Datastr.Resample.Sych.Average.InsoleGRFData.ave_l(:, 5)',...
     Datastr.Resample.Sych.Average.InsoleGRFData.std_l(:, 5)', 'b-', 0.6);
hold off
title('L CoPx')
ylabel('m')

subplot(2, 3, 3)
errorfill(phase,...
     Datastr.Resample.Sych.Average.ForcePlateGRFDataInCalcn.ave_l(:, 7)',...
     Datastr.Resample.Sych.Average.ForcePlateGRFDataInCalcn.std_l(:, 7)', 'r-', 1.0);
hold on
errorfill(phase,...
     Datastr.Resample.Sych.Average.InsoleGRFData.ave_l(:, 7)',...
     Datastr.Resample.Sych.Average.InsoleGRFData.std_l(:, 7)', 'b-', 0.6);
legend('FP-std', 'FP-mean', 'Insole-std', 'Insole-mean')
title('L CoPz')
ylabel('m')

subplot(2, 3, 4)
errorfill(phase,...
     Datastr.Resample.Sych.Average.ForcePlateGRFDataInCalcn.ave_r(:, 9)',...
     Datastr.Resample.Sych.Average.ForcePlateGRFDataInCalcn.std_r(:, 9)', 'r-', 1.0);
hold on;
errorfill(phase,...
     Datastr.Resample.Sych.Average.InsoleGRFData.ave_r(:, 9)',...
     Datastr.Resample.Sych.Average.InsoleGRFData.std_r(:, 9)', 'b-', 0.6);
hold off
title('R Fy')
ylabel('N')
xlabel('Phase (%)')

subplot(2, 3, 5)
errorfill(phase,...
     Datastr.Resample.Sych.Average.ForcePlateGRFDataInCalcn.ave_r(:, 11)',...
     Datastr.Resample.Sych.Average.ForcePlateGRFDataInCalcn.std_r(:, 11)', 'r-', 1.0);
hold on;
errorfill(phase,...
     Datastr.Resample.Sych.Average.InsoleGRFData.ave_r(:, 11)',...
     Datastr.Resample.Sych.Average.InsoleGRFData.std_r(:, 11)', 'b-', 0.6);
hold off
title('R CoPx')
ylabel('m')
xlabel('Phase (%)')

subplot(2, 3, 6)
errorfill(phase,...
     Datastr.Resample.Sych.Average.ForcePlateGRFDataInCalcn.ave_r(:, 13)',...
     Datastr.Resample.Sych.Average.ForcePlateGRFDataInCalcn.std_r(:, 13)', 'r-', 1.0);
hold on;
errorfill(phase,...
     Datastr.Resample.Sych.Average.InsoleGRFData.ave_r(:, 13)',...
     Datastr.Resample.Sych.Average.InsoleGRFData.std_r(:, 13)', 'b-', 0.6);
hold off
title('R CoPz')
ylabel('m')
xlabel('Phase (%)')
sgtitle(strcat('GRF - ', Datastr.Info.Trial))
savefig(strcat(Datastr.Info.SubjRoot, '\', fig_save_path, '\GRF-', Datastr.Info.Trial, '.fig'))

EMG_id = ["Sol","MGas","LGas","TA","Semi","BiFemoris","VasL","ReFemoris","VasM"];
figure()   % plot emgs
for e = 1:9
    subplot(5, 2, e)
     errorfill(phase,...
     Datastr.Resample.Sych.Average.EMG.ave_r(:, e+1)',...
     Datastr.Resample.Sych.Average.EMG.std_r(:, e+1)', 'r-', 1.0);
    title(EMG_id(e))
    ylabel('% MVC')
    if e == 9 || e == 10
        xlabel('Phase (%)')
    end
    if e == 2
        legend('EMG')
    end
    xlim([0, 101])
end

sgtitle(strcat('Muscle Excitation - ', Datastr.Info.Trial))
savefig(strcat(Datastr.Info.SubjRoot, '\', fig_save_path, '\MuscleExcitation-', Datastr.Info.Trial, '.fig'))
    
end