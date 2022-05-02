function hsMatrix = getHS(Fy, motionType)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the heel strike index matrix of the gait cyle from vertical force
    %
    if strcmpi(motionType, 'walk')
        threshold = 50;
    else
        threshold = 100;
    end
    index = zeros(length(Fy), 1);

    index(Fy > threshold) = 1;

    heelStrikes = find(diff(index, 1) == 1) - 2;  % two frames were compensated for the heel-strike
    
    % if the first strike is at the first 2 points, it will become negative
    % due to the compensation, then assign 1 to it.
    if heelStrikes(1) < 1   
        heelStrikes(1) = 1;
    end
    
    
    numStrike = length(heelStrikes);
    hsMatrix = [heelStrikes(1:end-1),  heelStrikes(2:end)];

    goodFlag = 0;
%     goodFy = 0;
        
    fprintf('Click the heel strike duration points that are not reasonable for the average gait cycle.\n')
    fprintf('If you think all the points are reasonable, just click the green dot in the middle left;\n')
    fprintf('If there are uncontinuely points too low or too high, then click the middle white part of the figure with your vertical line close to the point\n')
    fprintf('If there are continuely points too large starting from one point, click the top grey part of the figure with your vertical line colse to the first abnormal point;\n')
    fprintf('If there are continuely points too low starting from one point, then click the bottom grey part of the figure with your vertical line close to the first abnormal point;\n')
    fprintf('Only one click is acceptted in each iteration, the figure will be updated after your click, you can continuously working on the updated plot, or click the gree cirlcle if you want to finish\n');

    while ~goodFlag % && ~goodFy
        
       vRange = [min(hsMatrix(:, 2) - hsMatrix(:, 1)), max(hsMatrix(:, 2) - hsMatrix(:, 1))];
       
       fig1 = figure();
       subplot(2,1,1)
       rectangle('Position',[1,vRange(1),numStrike,(vRange(2)-vRange(1))/3],'FaceColor',[200 200 200]/255, 'LineWidth',1)
       hold on
       rectangle('Position',[1,vRange(2)-(vRange(2)-vRange(1))/3,numStrike,(vRange(2)-vRange(1))/3],'FaceColor',[200 200 200]/255, 'LineWidth',1)   
       hold on
       plot(hsMatrix(:, 2) - hsMatrix(:, 1), 'o'); 
       hold on
       plot(0, mean(vRange), 'go')
       title('Refer to the explination in the commend window for taking actions')
       hold off
       subplot(2,1,2)
       for i = 1:size(hsMatrix, 1)
           plot(1001:1001 + hsMatrix(i, 2) - hsMatrix(i, 1), Fy(hsMatrix(i, 1):hsMatrix(i, 2)));
           hold on
       end
       hold off
       
       [x,y]=ginput(1);
       
       FyError = 1e4;
       
       if x > 1000  % remove the gait cycle directly by clicking the Fy line
           for i = 1:size(hsMatrix, 1)
                FyInterp = Fy(hsMatrix(i, 1) + floor(x)-1001)*(1-x+floor(x)) +...
                           Fy(hsMatrix(i, 1) + ceil(x)-1001)*(1+x-ceil(x));
                
                if abs(FyInterp - y) <= FyError
                    FyError = abs(FyInterp - y);
                    Error_trial = i;
                end
           end
           hsMatrix(Error_trial, :) = [];
           numStrike = numStrike - 1;
       else  % remove the gait cycle directly by clicking the gati duration

           if sqrt(x^2 + (y-mean(vRange))^2) < 0.5
               goodFlag = 1;

           elseif y > vRange(2)-(vRange(2)-vRange(1))/3
               hsMatrix1 = hsMatrix(:, 1);
               hsMatrix2 = hsMatrix(:, 2);
               hsMatrix1(round(x)) = [];
               hsMatrix2(end) = [];
               hsMatrix = [hsMatrix1, hsMatrix2];
               numStrike = numStrike - 1;
           elseif y < vRange(2)-2*(vRange(2)-vRange(1))/3
               hsMatrix1 = hsMatrix(:, 1);
               hsMatrix2 = hsMatrix(:, 2);
               hsMatrix1(end) = [];
               hsMatrix2(round(x)) = [];
               hsMatrix = [hsMatrix1, hsMatrix2];
               numStrike = numStrike - 1;
           else
               hsMatrix(round(x), :) = [];
               numStrike = numStrike - 1;
           end
           
       end
       close(fig1)
    end
    
end

