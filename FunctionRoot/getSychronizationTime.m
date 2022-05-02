function data2_df = getSychronizationTime(data1, data2, TrialName)
%   This function calculate the sychronization time between the two dataset
%   motions, using the least square method
%   INPUTS:
%       data1: the first dataset M x N1 matrix
%       data2: the second dataset M x N2 matrix
%   

%% first make sure the sign are the same for two systems

[r1, c1] = size(data1);
[r2, c2] = size(data2);

if r1 < c1  % if row number is less than the column number, then rotate
    data1 = data1';
    [r1, c1] = size(data1);
end

if r2 < c2  % if row number is less than the column number, then rotate
    data2 = data2';
    [r2, c2] = size(data2);
end

if c1 ~= c2
    error('data1 and data2 should have the same number of variables')
end

rmin = min(r1, r2);
t = 1:size(data2(100:rmin-100, :), 1);

r = 0; % calculate the r
df_init = 0;

r_threshold = 0.7;
if strcmp(TrialName, 'lunge')
    r_threshold = 0.5;
end

while r < r_threshold % calculate the r
    % then apply the least square optimization

    data2_df = fmincon(@(df)calLsFit(df, data1(100:rmin-100, :),...
               data2(100:rmin-100, :), t), df_init, [], [], [], [], 0, 50);
    
    data2_sych = interp1(t, data2(100:rmin-100, :), t+data2_df, 'linear', 'extrap');

    r = trace(corr(data1(100:rmin-100, :), data2_sych))/c1; % calculate the r
    
    if df_init == 50
        error('Not able to get a fit of r > 0.7')
    end
    
    df_init = df_init + 5;

end
           

    function lsfit = calLsFit(df, data_tar, data_fit, t)  % the least square fit between x and y

        y_sych = interp1(t, data_fit, t+df, 'linear', 'extrap');

        lsfit = sum(sum((data_tar - y_sych).^2))/size(data_tar, 1) + 0.1*df^2;
        
    end

end

