function [Ang] = R2EA(R)
% Rotation matrix to Euler angles
%
% Input : Rotation matrix
% Output : Matrix of Euler angles in radian 
%          Each column is a solution, each row a rotation angle (alpha,beta,gamma)
%
% Assume given a rotation matrix R of size [3x3]
% Assume alpha, beta and gamma are the angles rotated about x,y and z axes
% Assume R = Rz*Ry*Rx ( = first rotation is about x axis, then y, then z),
% then:
% 
% R =
% [ cos(b)*cos(c), cos(c)*sin(a)*sin(b) - cos(a)*sin(c), sin(a)*sin(c) + cos(a)*cos(c)*sin(b)]
% [ cos(b)*sin(c), cos(a)*cos(c) + sin(a)*sin(b)*sin(c), cos(a)*sin(b)*sin(c) - cos(c)*sin(a)]
% [       -sin(b),                        cos(b)*sin(a),                        cos(a)*cos(b)]
%
%
% If -sin(b) ~= 1 or -1 then:
% you can solve b from R(3,1)
% you can solve a from R(3,2) and R(3,3), since (sin(a)*cos(b))/(cos(a)*cos(b)) = tan(a)
% you can solve c from R(2,1) and R(1,1), since (cos(b)*sin(c))/(cos(b)*cos(c)) = tan(c)
% Else
% These solutions do not hold, would lead to atan(0)
% But then R(1,3) = R(2,2) and -R(1,2) = R(2,3), and you can solve
%
% Mark Vlutters - Enschede - 2012

if abs(R(3,1)) ~= 1
    % Multiple solutions exist (depending which way you rotate)
    % Note that for atan2 the sign matters! 
    % e.g: To 'retrieve' the sign of sin(a) and cos(a) you need to divide by cos(b)
    
    beta_1 = asin(-R(3,1));
    beta_2 = pi - beta_1;
    
    alpha_1 = atan2(R(3,2)./cos(beta_1),R(3,3)./cos(beta_1));
    alpha_2 = atan2(R(3,2)./cos(beta_2),R(3,3)./cos(beta_2));
    
    gamma_1 = atan2(R(2,1)./cos(beta_1),R(1,1)./cos(beta_1));
    gamma_2 = atan2(R(2,1)./cos(beta_2),R(1,1)./cos(beta_2));
    
elseif R(3,1)==-1
    gamma_1 = 0; % actually, you can't solve gamma but that doesn't matter since alpha and beta provide sufficient info
    beta_1 = -pi/2;
    alpha_1 = gamma + atan(R(1,2),R(1,3));
elseif R(3,1)== 1
    gamma_1 = 0;
    beta_1 = pi/2;
    alpha_1 = -gamma + atan(-R(1,2),-R(1,3));
end

Ang = [ alpha_1 , alpha_2 ; beta_1 , beta_2 ; gamma_1 , gamma_2 ];

end

