function copConp = removeCoPSpeed(cop, speed, frame_rate)
% remove speed in the treadmill CoP X
        CoPx_LFP = cop(:, 1);
        CoPy_LFP = cop(:, 2);
        CoPx_RFP = cop(:, 3);
        CoPy_RFP = cop(:, 4);

        iframel = 0;
        iframer = 0;
        for i = 1:length(CoPx_LFP)-1
            if abs(CoPx_LFP(i+1) - CoPx_LFP(i)) > 1e-5
                iframel = iframel + 1;
                CoPx_LFP(i) = CoPx_LFP(i) + iframel/frame_rate*speed;
            else
                iframel = 0;
            end

            if abs(CoPx_RFP(i+1) - CoPx_RFP(i)) > 1e-5
                iframer = iframer + 1;
                CoPx_RFP(i) = CoPx_RFP(i) + iframer/frame_rate*speed;
            else
                iframer = 0;
            end            
        end

        ave_ly = mean(CoPy_LFP(CoPy_LFP~=0));
        ave_ry = mean(CoPy_RFP(CoPy_RFP~=0));

        CoPy_LFP = CoPy_LFP - ave_ly;
        CoPy_RFP = CoPy_RFP - ave_ry;
        
        copConp = [CoPx_LFP, CoPy_LFP, CoPx_RFP, CoPy_RFP];

end