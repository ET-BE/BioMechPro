function [interpVar_l, aveVar_l, stdVar_l, interpVar_r, aveVar_r, stdVar_r] = getGaitAve_noCheck(hsMatrix_l, hsMatrix_r, var, fRate, varName, figSavePath, plotGen)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the heel strike index matrix of the gait cyle from vertical force
   
    numStrikes_l = length(hsMatrix_l(:, 1));
    numStrikes_r = length(hsMatrix_r(:, 1));
    numVar = length(var(1, :));
    interpVar_l = zeros(100, numStrikes_l, numVar);
    interpVar_r = zeros(100, numStrikes_r, numVar);

    if strcmpi(varName, 'IKAngData') || strcmpi(varName, 'IMUAngData')
        
        for i = 1:numStrikes_l
            time_org = (hsMatrix_l(i, 1):hsMatrix_l(i, 2))/fRate;
            phases = linspace(hsMatrix_l(i, 1)/fRate, hsMatrix_l(i, 2)/fRate, 100);
            interpVar_l(:, i, :) = interp1(time_org, var(hsMatrix_l(i, 1):hsMatrix_l(i, 2), :), phases);
        end
        
        for i = 1:numStrikes_r
            time_org = (hsMatrix_r(i, 1):hsMatrix_r(i, 2))/fRate;
            phases = linspace(hsMatrix_r(i, 1)/fRate, hsMatrix_r(i, 2)/fRate, 100);
            interpVar_r(:, i, :) = interp1(time_org, var(hsMatrix_r(i, 1):hsMatrix_r(i, 2), :), phases);
        end
        
        joint_labels = ["LHip", "RHip", "LKnee", "RKnee", "LAnkle", "RAnkle"];
        joint_angle_id = [15, 8, 18, 11, 19, 12];
                
        if strcmpi(plotGen, 'True')
            fig1 = figure();   % plot joint angles
            for j = 1:6
                subplot(3, 2, j)
                if rem(j, 2)
                    plot(1:100,...
                    interpVar_l(:, :, joint_angle_id(j)), 'linewidth', 1);
                else
                    plot(1:100,...
                    interpVar_r(:, :, joint_angle_id(j)), 'linewidth', 1);
                end
                hold off;
                title(joint_labels(j))
                if j == 5 || j == 6
                    xlabel('Phase (%)')
                end
                if j == 1|| j == 3|| j == 5
                    ylabel('Angle (deg.)')
                end
            end
            sgtitle(strcat('Joint Angle - ', varName))
            savefig(strcat(figSavePath, '\JointAngle_', varName, '.fig'))
            close(fig1)
        end
        
    elseif strcmpi(varName, 'IDTrqData') || strcmpi(varName, 'IDTrqData_Portable')
        
        for i = 1:numStrikes_l
            time_org = (hsMatrix_l(i, 1):hsMatrix_l(i, 2))/fRate;
            phases = linspace(hsMatrix_l(i, 1)/fRate, hsMatrix_l(i, 2)/fRate, 100);
            interpVar_l(:, i, :) = interp1(time_org, var(hsMatrix_l(i, 1):hsMatrix_l(i, 2), :), phases);
        end
        
        for i = 1:numStrikes_r
            time_org = (hsMatrix_r(i, 1):hsMatrix_r(i, 2))/fRate;
            phases = linspace(hsMatrix_r(i, 1)/fRate, hsMatrix_r(i, 2)/fRate, 100);
            interpVar_r(:, i, :) = interp1(time_org, var(hsMatrix_r(i, 1):hsMatrix_r(i, 2), :), phases);
        end
                
        joint_labels = ["LHip", "RHip", "LKnee", "RKnee", "LAnkle", "RAnkle"];
        joint_torque_id = [11, 8, 18, 17, 20, 19];
        
        if strcmpi(plotGen, 'True')
            fig1 = figure();  % plot joint torques
            for j = 1:6
                subplot(3, 2, j)
                if rem(j, 2)
                    plot(1:100,...
                    interpVar_l(:, :, joint_torque_id(j)), 'linewidth', 1);
                else
                    plot(1:100,...
                    interpVar_r(:, :, joint_torque_id(j)), 'linewidth', 1);
                end
                hold off
                title(joint_labels(j))
                if j == 5 || j == 6
                    xlabel('Time (s)')
                end
                if j == 1|| j == 3|| j == 5
                    ylabel('Torque (Nm)')
                end
            end
            sgtitle(strcat('Joint Torque - ', varName))
            savefig(strcat(figSavePath, '\JointTorque_', varName, '.fig'))
            close(fig1)
        end
        
    elseif strcmpi(varName, 'ForcePlateGRFData')
        
%         var(:, [5, 7, 11 ,13]) = removeCoPSpeed(var(:, [5, 7, 11 ,13]), speed, fRate);
        
        for i = 1:numStrikes_l
            time_org = (hsMatrix_l(i, 1):hsMatrix_l(i, 2))/fRate;
            phases = linspace(hsMatrix_l(i, 1)/fRate, hsMatrix_l(i, 2)/fRate, 100);
            interpVar_l(:, i, :) = interp1(time_org, var(hsMatrix_l(i, 1):hsMatrix_l(i, 2), :), phases);
            interpVar_l(:, i, 5:7) = interp1(time_org, var(hsMatrix_l(i, 1):hsMatrix_l(i, 2), 5:7) - var(hsMatrix_l(i, 1), 5:7), phases);
        end
        
        for i = 1:numStrikes_r
            time_org = (hsMatrix_r(i, 1):hsMatrix_r(i, 2))/fRate;
            phases = linspace(hsMatrix_r(i, 1)/fRate, hsMatrix_r(i, 2)/fRate, 100);
            interpVar_r(:, i, :) = interp1(time_org, var(hsMatrix_r(i, 1):hsMatrix_r(i, 2), :), phases);
            interpVar_r(:, i, 11:13) = interp1(time_org, var(hsMatrix_r(i, 1):hsMatrix_r(i, 2), 11:13) - var(hsMatrix_r(i, 1), 11:13), phases);
        end
       
        if strcmpi(plotGen, 'True')
            fig1 = figure();  % plot Fy and CoPs
            subplot(2, 3, 1)
            plot(1:100, interpVar_l(:, :, 3), 'linewidth', 1);
            title('L Fy')
            ylabel('N')

            subplot(2, 3, 2)
            plot(1:100, interpVar_l(:, :, 5), '.', 'linewidth', 1);
            title('L CoPx')
            ylabel('m')

            subplot(2, 3, 3)
            plot(1:100, interpVar_l(:, :, 7), '.', 'linewidth', 1);
            title('L CoPz')
            ylabel('m')

            subplot(2, 3, 4)
            plot(1:100, interpVar_r(:, :, 9), 'linewidth', 1);
            title('R Fy')
            ylabel('N')
            xlabel('Time (s)')

            subplot(2, 3, 5)
            plot(1:100, interpVar_r(:, :, 11), '.', 'linewidth', 1);
            title('R CoPx')
            ylabel('m')
            xlabel('Time (s)')

            subplot(2, 3, 6)
            plot(1:100, interpVar_r(:, :, 13), '.', 'linewidth', 1);
            title('R CoPz')
            ylabel('m')
            xlabel('Time (s)')
            sgtitle(strcat('GRF - ', varName))
            savefig(strcat(figSavePath, '\GRF_', varName, '.fig'))
            close(fig1)
        end

    elseif strcmpi(varName, 'InsoleGRFData')
        for i = 1:numStrikes_l
            time_org = (hsMatrix_l(i, 1):hsMatrix_l(i, 2))/fRate;
            phases = linspace(hsMatrix_l(i, 1)/fRate, hsMatrix_l(i, 2)/fRate, 100);
            interpVar_l(:, i, :) = interp1(time_org, var(hsMatrix_l(i, 1):hsMatrix_l(i, 2), :), phases);
        end
        
        for i = 1:numStrikes_r
            time_org = (hsMatrix_r(i, 1):hsMatrix_r(i, 2))/fRate;
            phases = linspace(hsMatrix_r(i, 1)/fRate, hsMatrix_r(i, 2)/fRate, 100);
            interpVar_r(:, i, :) = interp1(time_org, var(hsMatrix_r(i, 1):hsMatrix_r(i, 2), :), phases);
        end
        if strcmpi(plotGen, 'True')
            fig1 = figure();  % plot Fy and CoPs
            subplot(2, 3, 1)
            plot(1:100, interpVar_l(:, :, 3), 'linewidth', 1);
            title('L Fy')
            ylabel('N')

            subplot(2, 3, 2)
            plot(1:100, interpVar_l(:, :, 5), '.', 'linewidth', 1);
            title('L CoPx')
            ylabel('m')

            subplot(2, 3, 3)
            plot(1:100, interpVar_l(:, :, 7), '.', 'linewidth', 1);
            title('L CoPz')
            ylabel('m')

            subplot(2, 3, 4)
            plot(1:100, interpVar_r(:, :, 9), 'linewidth', 1);
            title('R Fy')
            ylabel('N')
            xlabel('Time (s)')

            subplot(2, 3, 5)
            plot(1:100, interpVar_r(:, :, 11), '.', 'linewidth', 1);
            title('R CoPx')
            ylabel('m')
            xlabel('Time (s)')

            subplot(2, 3, 6)
            plot(1:100, interpVar_r(:, :, 13), '.', 'linewidth', 1);
            title('R CoPz')
            ylabel('m')
            xlabel('Time (s)')
            sgtitle(strcat('GRF - ', varName))
            savefig(strcat(figSavePath, '\GRF_', varName, '.fig'))
            close(fig1)
        end
        
    elseif strcmpi(varName, 'ForcePlateGRFDataInCalcn')
                
        for i = 1:numStrikes_l
            time_org = (hsMatrix_l(i, 1):hsMatrix_l(i, 2))/fRate;
            phases = linspace(hsMatrix_l(i, 1)/fRate, hsMatrix_l(i, 2)/fRate, 100);
            interpVar_l(:, i, :) = interp1(time_org, var(hsMatrix_l(i, 1):hsMatrix_l(i, 2), :), phases);
        end
        
        for i = 1:numStrikes_r
            time_org = (hsMatrix_r(i, 1):hsMatrix_r(i, 2))/fRate;
            phases = linspace(hsMatrix_r(i, 1)/fRate, hsMatrix_r(i, 2)/fRate, 100);
            interpVar_r(:, i, :) = interp1(time_org, var(hsMatrix_r(i, 1):hsMatrix_r(i, 2), :), phases);
        end
       
        if strcmpi(plotGen, 'True')
            fig1 = figure();  % plot Fy and CoPs
            subplot(2, 3, 1)
            plot(1:100, interpVar_l(:, :, 3), 'linewidth', 1);
            title('L Fy')
            ylabel('N')

            subplot(2, 3, 2)
            plot(1:100, interpVar_l(:, :, 5), '.', 'linewidth', 1);
            title('L CoPx')
            ylabel('m')

            subplot(2, 3, 3)
            plot(1:100, interpVar_l(:, :, 7), '.', 'linewidth', 1);
            title('L CoPz')
            ylabel('m')

            subplot(2, 3, 4)
            plot(1:100, interpVar_r(:, :, 9), 'linewidth', 1);
            title('R Fy')
            ylabel('N')
            xlabel('Time (s)')

            subplot(2, 3, 5)
            plot(1:100, interpVar_r(:, :, 11), '.', 'linewidth', 1);
            title('R CoPx')
            ylabel('m')
            xlabel('Time (s)')

            subplot(2, 3, 6)
            plot(1:100, interpVar_r(:, :, 13), '.', 'linewidth', 1);
            title('R CoPz')
            ylabel('m')
            xlabel('Time (s)')
            sgtitle(strcat('GRF - ', varName))
            savefig(strcat(figSavePath, '\GRFInCalcn_', varName, '.fig'))
            close(fig1)
        end

    elseif strcmpi(varName, 'EMG')
        
        for i = 1:numStrikes_r
            time_org = (hsMatrix_r(i, 1):hsMatrix_r(i, 2))/fRate;
            phases = linspace(hsMatrix_r(i, 1)/fRate, hsMatrix_r(i, 2)/fRate, 100);
            interpVar_r(:, i, :) = interp1(time_org, var(hsMatrix_r(i, 1):hsMatrix_r(i, 2), :), phases);
        end
        
        if strcmpi(plotGen, 'True')
            EMG_id = ["Sol","MGas","LGas","TA","Semi","BiFemoris","VasL","ReFemoris","VasM"];
            fig1 = figure();
            for e = 1:9
                subplot(5, 2, e)
                plot(1:100,interpVar_r(:, :, e), 'linewidth', 0.5)
                title(EMG_id(e))
                ylabel('% MVC')
                if e == 9 || e == 10
                    xlabel('Time (s)')
                end
                if e == 2
                    legend('EMG')
                end
            end

            sgtitle(strcat('Muscle Excitation - ', varName))
            savefig(strcat(figSavePath, '\EMG_', varName, '.fig'))
            close(fig1)
        end
    end
    
    if strcmpi(varName, 'EMG')
    aveVar_l = squeeze(mean(interpVar_l, 2));
    aveVar_l = [(1:100)', aveVar_l];
    stdVar_l = squeeze(std(interpVar_l, 0, 2));
    stdVar_l = [zeros(100, 1), stdVar_l];
    
    aveVar_r = squeeze(mean(interpVar_r, 2));
    aveVar_r = [(1:100)', aveVar_r];
    stdVar_r = squeeze(std(interpVar_r, 0, 2));
    stdVar_r = [zeros(100, 1), stdVar_r];
    else
    aveVar_l = squeeze(mean(interpVar_l, 2));
    aveVar_l(:, 1) = 1:100;
    stdVar_l = squeeze(std(interpVar_l, 0, 2));
    stdVar_l(:, 1) = zeros(100, 1);
    
    aveVar_r = squeeze(mean(interpVar_r, 2));
    aveVar_r(:, 1) = 1:100;
    stdVar_r = squeeze(std(interpVar_r, 0, 2));
    stdVar_r(:, 1) = zeros(100, 1);
    end
    

end