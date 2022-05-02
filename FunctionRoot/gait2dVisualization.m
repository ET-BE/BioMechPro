
osim_joint = [2.129 ...
             -13.0 -0.43 38.807 2.18 0.939 0.351 ...
             9.391 1.827 1.092 ...
             -3.022 ...
             6.894 -1.962 0.0 ...
             7.218 -0.094 5.582 ...
             -1.005 ...
             6.266 1.306 0.0 ...
             16.686 3.887 -5.447];

osim_grf = [0 44 0 0 0 0 0 528 0 -0.15313 0 -0.09361];
com = [0.15 1.015 0.402];
height = 1.78;
lcolor = 'b';
forceplate = 0;
insole = 1;

figure()
haxes = axes;
gait2dVisualization1(osim_joint, osim_grf, com, height, haxes, lcolor, forceplate, insole)

function gait2dVisualization1(osim_joint, osim_grf, com, height, haxes, lcolor, forceplate, insole)
% This function is to generate the stick plot of the 2d gait model

% load plotting information
timeframe = osim_joint(1);
trunk_bending = osim_joint(end-2)*pi/180;
pelvis_x = osim_joint(5);
pelvis_y = osim_joint(6);
pelvis_tilt = osim_joint(2)*pi/180;
rhip_flex = osim_joint(8)*pi/180;
lhip_flex = osim_joint(15)*pi/180;
rknee_flex = osim_joint(11)*pi/180;
lknee_flex = osim_joint(18)*pi/180;
rfoot_flex = osim_joint(12)*pi/180;
lfoot_flex = osim_joint(19)*pi/180;

lFy = osim_grf(2);
lPx = osim_grf(4);
rFy = osim_grf(8);
rPx = osim_grf(10);

comx = com(1);
comy = com(2);

% set body demensions based on the body height
ltrunk = height*0.288;
lthigh = height*0.245;
lshank = height*0.246;
hfoot  = height*0.039;
lfoot  = height*0.152;
ht_ra = 0.3;


% correct the pelvis rotation
rhip_flex = rhip_flex + pelvis_tilt;
lhip_flex = lhip_flex + pelvis_tilt;

% generate the joint points list based on the joint angles
pointx = [pelvis_x - ltrunk*sin(trunk_bending),...  % neck point
          pelvis_x, pelvis_x + lthigh*sin(rhip_flex),...  % pelvis location & rknee point
          pelvis_x + lthigh*sin(rhip_flex) + lshank*sin(rknee_flex + rhip_flex),...  % rankle point
          pelvis_x + lthigh*sin(rhip_flex) + lshank*sin(rknee_flex ...  
          + rhip_flex)+ hfoot*sin(rknee_flex + rhip_flex + rfoot_flex)...  % rheel point
          - ht_ra*lfoot*cos(rknee_flex + rhip_flex + rfoot_flex),...
          pelvis_x + lthigh*sin(rhip_flex) + lshank*sin(rknee_flex ...
          + rhip_flex)+ hfoot*sin(rknee_flex + rhip_flex + rfoot_flex)...  % rtoe point
          + (1-ht_ra)*lfoot*cos(rknee_flex + rhip_flex + rfoot_flex),...
          pelvis_x + lthigh*sin(rhip_flex) + lshank*sin(rknee_flex + rhip_flex),... % rankle point
          pelvis_x + lthigh*sin(lhip_flex),...  % lknee point
          pelvis_x + lthigh*sin(lhip_flex) + lshank*sin(lknee_flex + lhip_flex),...  % lankle point
          pelvis_x + lthigh*sin(lhip_flex) + lshank*sin(lknee_flex ...  
          + lhip_flex)+ hfoot*sin(lknee_flex + lhip_flex + lfoot_flex)...  % lheel point
          - ht_ra*lfoot*cos(lknee_flex + lhip_flex + lfoot_flex),...
          pelvis_x + lthigh*sin(lhip_flex) + lshank*sin(lknee_flex ...
          + lhip_flex)+ hfoot*sin(lknee_flex + lhip_flex + lfoot_flex)...  % ltoe point
          + (1-ht_ra)*lfoot*cos(lknee_flex + lhip_flex + lfoot_flex),...
          pelvis_x + lthigh*sin(lhip_flex) + lshank*sin(lknee_flex + lhip_flex)... % lankle point
          ];
pointy = [pelvis_y + ltrunk*cos(trunk_bending),...  % neck point
          pelvis_y, pelvis_y - lthigh*cos(rhip_flex),...         % pelvis location & rknee point
          pelvis_y - lthigh*cos(rhip_flex) - lshank*cos(rknee_flex + rhip_flex),...  % rankle point
          pelvis_y - lthigh*cos(rhip_flex) - lshank*cos(rknee_flex ...  
          + rhip_flex)- hfoot*cos(rknee_flex + rhip_flex + rfoot_flex)...  % rheel point
          - ht_ra*lfoot*sin(rknee_flex + rhip_flex + rfoot_flex),...
          pelvis_y - lthigh*cos(rhip_flex) - lshank*cos(rknee_flex ...
          + rhip_flex)- hfoot*cos(rknee_flex + rhip_flex + rfoot_flex)...  % rtoe point
          + (1-ht_ra)*lfoot*sin(rknee_flex + rhip_flex + rfoot_flex),...
          pelvis_y - lthigh*cos(rhip_flex) - lshank*cos(rknee_flex + rhip_flex),... % rankle point
          pelvis_y - lthigh*cos(lhip_flex),...  % lknee point
          pelvis_y - lthigh*cos(lhip_flex) - lshank*cos(lknee_flex + lhip_flex),...  % lankle point
          pelvis_y - lthigh*cos(lhip_flex) - lshank*cos(lknee_flex ...  
          + lhip_flex) - hfoot*cos(lknee_flex + lhip_flex + lfoot_flex)...  % lheel point
          - ht_ra*lfoot*sin(lknee_flex + lhip_flex + lfoot_flex),...
          pelvis_y - lthigh*cos(lhip_flex) - lshank*cos(lknee_flex ...
          + lhip_flex) - hfoot*cos(lknee_flex + lhip_flex + lfoot_flex)...  % ltoe point
          + (1-ht_ra)*lfoot*sin(lknee_flex + lhip_flex + lfoot_flex),...
          pelvis_y - lthigh*cos(lhip_flex) - lshank*cos(lknee_flex + lhip_flex)... % lankle point
          ];

% generate the Fy force vector
if insole
    rfx = pelvis_x + lthigh*sin(rhip_flex) + lshank*sin(rknee_flex ...  
        + rhip_flex)+ hfoot*sin(rknee_flex + rhip_flex + rfoot_flex)...  % rCoPx points
        + (ht_ra + rPx)*lfoot*cos(rknee_flex + rhip_flex + rfoot_flex);
    
    rfy = pelvis_y - lthigh*cos(rhip_flex) - lshank*cos(rknee_flex ...
          + rhip_flex)- hfoot*cos(rknee_flex + rhip_flex + rfoot_flex)...  % rtoe point
          + (ht_ra + rPx)*lfoot*sin(rknee_flex + rhip_flex + rfoot_flex);
      
    lfx = pelvis_x + lthigh*sin(lhip_flex) + lshank*sin(lknee_flex ...  
        + lhip_flex)+ hfoot*sin(lknee_flex + lhip_flex + lfoot_flex)...  % rCoPx points
        + (ht_ra + lPx)*lfoot*cos(lknee_flex + lhip_flex + lfoot_flex);
    
    lfy = pelvis_y - lthigh*cos(lhip_flex) - lshank*cos(lknee_flex ...
          + lhip_flex)- hfoot*cos(lknee_flex + lhip_flex + lfoot_flex)...  % rtoe point
          + (ht_ra + lPx)*lfoot*sin(lknee_flex + lhip_flex + lfoot_flex);
    
    fx = [rfx, rfx, lfx, lfx];
    fy = [rfy, rfy + rFy/1000, lfy, lfy + lFy/1000];
    
elseif forceplate
    fx = [rPx, rPx, lPx, lPx];
    fy = [0, rFy/1000, 0, lFy/1000];
    
else
    error('Force type cannot be both force plate and insole')
end
      

% generate the plotting lines
lw = 2;
ms = 6;
plot(haxes, pointx(1:7), pointy(1:7), 'o-', 'linewidth', lw, 'MarkerSize', ms, 'Color', lcolor);
hold on
plot(haxes, pointx([2, 8]), pointy([2, 8]), 'o-', 'linewidth', lw,  'MarkerSize', ms, 'Color', lcolor);
hold on
plot(haxes, pointx(8:end), pointy(8:end), 'o-', 'linewidth', lw,  'MarkerSize', ms, 'Color', lcolor);
hold on;
plot(haxes, comx, comy, 'o',  'MarkerSize', ms+6, 'Color', lcolor)
hold on;
plot(haxes, fx(1:2), fy(1:2), '-', 'linewidth', lw+1, 'Color', lcolor)
hold on
plot(haxes, fx(3:4), fy(3:4), '-', 'linewidth', lw+1, 'Color', lcolor)

xl = xlim;
yl = ylim;
text(xl(1), yl(2), num2str(timeframe))

end