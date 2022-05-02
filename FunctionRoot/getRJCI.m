function [C3Ddata] = getRJCI(C3Ddata,sex,mass)
%% getRJCI (Rotations, Joints, COMs and Inertia)
% Get Rotation matrices and translations that transforms segment coordinates to global coordinates
% Get Joint positions
% Get COM positions per segment, and total COM
% Get Inertia properties per segment
% Get Mass properties per segment
% 
% See Dumas et al, 2006, Adjustments to McConville...
% 
% INPUT)
% C3Ddata structure, with at least the fields:
% 
% OUTPUT)
% C3Ddata structure, with added fields:
% 
% NOTES)
% All COM location percentages are rounded to 0 if it's below 1% in table 2.
% 
% The arms are not yet taken up in this file
% 
% seg2GlobD and jointData basically contain the same info, but they are 
% separate variables as you can find more joint centers (in jointData) than
% you have segments (stored in Seg2Glob). For example you can find the CJC
% and both HJC from the pelvis.
% 
% Note that the inertia needs a double transformation with seg2globR
% This is because the inertia arises from a quadratic distance
% Iglob = m .* Dglob*Dglob'
% If :
% Dglob = R*Dloc 
% Then :
% Iglob = m .* (R*Dloc * Dloc'*R')
% Iglob = R Iloc R'
% 
% The inertia tensor about the COM must be obtained by applying the
% parallel axis theorem to the inertia tensor.
% 
% Update history)
% 14-11-2014 : Mark Vlutters


%% Get information

% Predefined names of landmarks
ProbeCell = ...
    {'CAL','CM1L','CM5L', 'CAR','CM1R','CM5R', ...          % Feet
    'MLL','MML','CFL', 'MLR','MMR','CFR' ...                % Lower legs
    'ELL','EML', 'ELR','EMR', ...                           % Upper legs
    'SIASL','SIPSL', 'SIASR','SIPSR',...                    % Pelvis
    'C7','PX','IJ',...                                      % Trunk
    'OCC','HV','SEL'};                                      % Head

% Initiate variables
labelall = {};
labelallC = {};
markerall = [];

% Get all labels and data
if isfield(C3Ddata.Marker,'MarkerDataLabel')
    labelall(end+1:end+length(C3Ddata.Marker.MarkerDataLabel)) = C3Ddata.Marker.MarkerDataLabel;
    labelallC(end+1:end+length(C3Ddata.Marker.MarkerDataLabel)) = cell([1 length(C3Ddata.Marker.MarkerDataLabel)]);
    markerall(:,end+1:end+length(C3Ddata.Marker.MarkerDataLabel),:) = C3Ddata.Marker.MarkerData;
end
if isfield(C3Ddata.Marker,'ProbedDataLabel');
    labelall(end+1:end+length(C3Ddata.Marker.ProbedDataLabel)) = C3Ddata.Marker.ProbedDataLabel;
    labelallC(end+1:end+length(C3Ddata.Marker.ProbedDataLabel)) = C3Ddata.Marker.ProbedDataLabelC;
    markerall(:,end+1:end+length(C3Ddata.Marker.ProbedDataLabel),:) = C3Ddata.Marker.ProbedData;  % NOTE: markerall is not unused, it's used in an EVAL statement
end


% Compare names in ProbeCell with names in the available label cells.
% If a marker with that name (in it) is there, create a variable with the name from ProbeCell
for ilabel = 1:length(labelall)
    % regexpi is used here instead of strcmpi, because the names might not exactly match with those in ProbeCell (e.g. in case of preceding with cluster name)
    if any( cellfun('isempty' , regexpi(labelall(ilabel),ProbeCell)) == 0 ) 
        
        % Get the right name for the variable
        varname = ProbeCell{ cellfun('isempty' , regexpi(labelall(ilabel),ProbeCell)) == 0 };

        % Creat the variable (note that the name may NOT occur multiple times, i.e. varname may only have one value True)
        eval( [varname ' = squeeze( markerall(:,' num2str(ilabel) ',:) );'] );
        
    end
end


%% See what is possible with the acquired data
vari = who;

isPELVIS = sum(ismember({'SIASL','SIASR','SIPSL','SIPSR'},vari)) == 4;
isTRUNK = sum(ismember({'SIASL','SIASR','SIPSL','SIPSR','C7','IJ','PX'},vari)) == 7;  % NOTE: Shoulder joint centers are not estimated here !
isHEAD = sum(ismember({'SIASL','SIASR','SIPSL','SIPSR','C7','IJ','PX','SEL','HV'},vari)) == 9;

isLUPPER = sum(ismember({'SIASL','SIASR','SIPSL','SIPSR','ELL','EML'},vari)) == 6;
isRUPPER = sum(ismember({'SIASL','SIASR','SIPSL','SIPSR','ELR','EMR'},vari)) == 6;

isLLOWER = sum(ismember({'MLL','MML','CFL','ELL','EML'},vari)) == 5;
isRLOWER = sum(ismember({'MLR','MMR','CFR','ELR','EMR'},vari)) == 5;

isLFOOT = sum(ismember({'MLL','MML','CAL','CM1L','CM5L'},vari)) == 5;
isRFOOT = sum(ismember({'MLR','MMR','CAR','CM1R','CM5R'},vari)) == 5;

isAJCL = isLFOOT || isLLOWER;
isAJCR = isRFOOT || isRLOWER;
isKJCL = isLLOWER || isLUPPER;
isKJCR = isRLOWER || isRUPPER;

% Number of segments
nseg = sum([isLFOOT isRFOOT isLLOWER isRLOWER isLUPPER isRUPPER isPELVIS isTRUNK isHEAD]);

% Number of joints
njoint = sum([isAJCL isAJCR isKJCL isKJCR 3.*isPELVIS isTRUNK]);

seg2globR = zeros( eval( ['size( ' varname ',1)'] ) , nseg , 3 , 3 );
seg2globD = zeros( eval( ['size( ' varname ',1)'] ) , nseg , 3 );
comData = zeros( eval( ['size( ' varname ',1)'] ) , nseg , 3 );
jointData = zeros( eval( ['size( ' varname ',1)'] ) , njoint , 3);
segInertia = zeros( eval( ['size( ' varname ',1)'] ) , nseg , 3 , 3 );
segInertiaO = zeros( eval( ['size( ' varname ',1)'] ) , nseg , 3 , 3 );

% Initiate variables
n = 0; % Segment counter
m = 0; % Joint counter


%% CENTER
%% PELVIS

if isPELVIS
    n = n + 1;
    
    % z-axis : from SIASL to SIASR
    zvec = SIASR - SIASL;
    zvecNorm = zvec ./ repmat( sqrt(sum(zvec.^2,2)) , [1 3] );

    % y-axis : normal to the plane containing SIASL, SIASR and the midpoint between SIPSL and SIPSR
    yvec = cross( SIASR - (SIPSR + SIPSL)./2 , SIASL - (SIPSR + SIPSL)./2 ); % Points cranial by using cross(SIASR,SIASL)
    yvecNorm = yvec ./ repmat( sqrt(sum(yvec.^2,2)) , [1 3] );
    
    % x-axis : cross-product of y and z axes
    xvecNorm = cross(yvecNorm,zvecNorm);
    
    % Get rotation matrix to get from a segment coordinates to global coordinates
    seg2globR(:,n,1:3,1) = xvecNorm;
    seg2globR(:,n,1:3,2) = yvecNorm;
    seg2globR(:,n,1:3,3) = zvecNorm;
        
    % Pelvis width: (note: should be approximately the same value for all samples)
    pelvWid = sqrt( zvec(:,1).^2 + zvec(:,2).^2 + zvec(:,3).^2 );
    
    % Origin (Cappozzo)
    orig = (SIASR + SIASL)./2;
    
    % Estimate HJC and LJC (see Dumas, Appendix B)
    HJCL = zeros( size(seg2globR,1) , 3 );
    HJCR = zeros( size(seg2globR,1) , 3 );
    LJC = zeros( size(seg2globR,1) , 3 );
    if strcmpi(sex,'m')
            % Mass
            segMass(n) = 0.142 .* mass;
            
            % COM local
            comLocal = [0.028 , -0.280 , 0]';
        
            % Inertia matrix 
            iMat = [1.01 -0.25 -0.12 ; -0.25 1.06 -0.08 ; -0.12 -0.08 0.95];
            
        for i = 1:size(seg2globR,1)
            % Joints
            HJCL(i,:) = orig(i,:)' + reshape(seg2globR(i,n,:,:),[3 3]) * [-0.208*pelvWid(i) , -0.278*pelvWid(i) , -0.361*pelvWid(i)]' ;
            HJCR(i,:) = orig(i,:)' + reshape(seg2globR(i,n,:,:),[3 3]) * [-0.208*pelvWid(i) , -0.278*pelvWid(i) ,  0.361*pelvWid(i)]' ;
            LJC(i,:)  = orig(i,:)' + reshape(seg2globR(i,n,:,:),[3 3]) * [-0.264*pelvWid(i) ,  0.126*pelvWid(i) , 0]' ;
            
            % Segment length : LJC to projection hip joints in sagittal plane (assume symmetric)
            segLen = sqrt( sum( ( LJC(i,:) - (HJCL(i,:) + HJCR(i,:))./2 ).^2 , 2) );
    
            % COM global
            comData(i,n,:) = LJC(i,:)' + reshape(seg2globR(i,n,:,:),[3 3]) * (segLen.*comLocal);
            
            % Inertia global
            segInertiaO(i,n,:,:) = reshape(seg2globR(i,n,:,:),[3 3]) * ( segMass(n).* (iMat*segLen).^2 ) * reshape(seg2globR(i,n,:,:),[3 3])'; % About origin
            parAxShft = eye(3).*norm(segLen.*comLocal).^2 - (segLen.*comLocal)*(segLen.*comLocal.');
            segInertia(i,n,:,:) = reshape(seg2globR(i,n,:,:),[3 3]) * ( segMass(n).* ((iMat*segLen).^2 + parAxShft) ) * reshape(seg2globR(i,n,:,:),[3 3])'; % About COM
           
        end
    elseif strcmpi(sex,'f')
            % Mass
            segMass(n) = 0.146 .* mass;
            
            % COM local
            comLocal = [0 -0.232 0]';
            
            % Inertia matrix 
            iMat = [0.91 -0.34 -0.01 ; -0.34 1 -0.01 ; -0.01 -0.01 0.79];
            
        for i = 1:size(seg2globR,1)
            % Joints
            HJCL(i,:) = orig(i,:)' + reshape(seg2globR(i,n,:,:),[3 3]) * [-0.197*pelvWid(i) , -0.27*pelvWid(i) , -0.372*pelvWid(i)]' ;
            HJCR(i,:) = orig(i,:)' + reshape(seg2globR(i,n,:,:),[3 3]) * [-0.197*pelvWid(i) , -0.27*pelvWid(i) ,  0.372*pelvWid(i)]' ;
            LJC(i,:)  = orig(i,:)' + reshape(seg2globR(i,n,:,:),[3 3]) * [-0.289*pelvWid(i) ,  0.172*pelvWid(i) , 0]' ;
            
            % Segment length : LJC to projection hip joints in sagittal plane (assume symmetric)
            segLen = sqrt( sum( ( LJC(i,:) - (HJCL(i,:) + HJCR(i,:))./2 ).^2 , 2) );
            
            % COM global
            comData(i,n,:) = LJC(i,:)' + reshape(seg2globR(i,n,:,:),[3 3]) * (segLen.*comLocal); % NOTE: LJC is the origin in Dumas, see below
            
            % Inertia global
            segInertiaO(i,n,:,:) = reshape(seg2globR(i,n,:,:),[3 3]) * ( segMass(n).* (iMat*segLen).^2 ) * reshape(seg2globR(i,n,:,:),[3 3])'; % About origin
            parAxShft = eye(3).*norm(segLen.*comLocal).^2 - (segLen.*comLocal)*(segLen.*comLocal.');
            segInertia(i,n,:,:) = reshape(seg2globR(i,n,:,:),[3 3]) * ( segMass(n).* ((iMat*segLen).^2 + parAxShft) ) * reshape(seg2globR(i,n,:,:),[3 3])'; % About COM
        end
    end
    
    % Origin (Dumas)
    % NOTE: what happens here is a bit odd, but Dumas uses Cappozzo's origin to estimate LJC and HJC, then uses the LJC as origin
    seg2globD(:,n,:) = LJC;
    
    % Store joints
    m = m + 1;
    jointData(:,m,:) = LJC;
    jointDataLabel{m} = 'LJC';
    m = m + 1;
    jointData(:,m,:) = HJCL;
    jointDataLabel{m} = 'HJCL';
    m = m + 1;
    jointData(:,m,:) = HJCR;
    jointDataLabel{m} = 'HJCR';
    
    % Store COM
    % See above: occurs inside loop
    
    % Obtain the names for COMLabelCell, using the cluster names
    comDataLabel{n} = labelallC{strcmpi('SIASL',labelall)};
    
end

%% TRUNK

if isTRUNK
    n = n + 1;
    
    % NOTE : requires the LJC from the PELVIS as origin
    
    % Get the CJC
    % As the CJC needs to be in the sagittal plane, you need a third point to define this plane, which is not mentioned in Dumas.
    % Here I use the PX
    % First I find the intersection of the vectors IJ-PX and the vector making an angle of 8 (14) deg with the vector C7-IJ
    % The CJC is along the vector Intersection-C7
    % Beta is the (smallest) angle between vectors C7-PX and IJ-PX
    % Gamma is the angle C7-Intersect-IJ
    if strcmpi(sex,'m')
        alpha = 8*pi/180;
        beta = acos( sum( (C7-IJ).*(IJ-PX),2 ) ./ ( sqrt( sum((C7-IJ).^2,2) ) .* sqrt( sum((IJ-PX).^2,2) ) ) );
        gamma = pi - alpha - beta;
        INTERSECT = IJ + repmat( ( sin(alpha)./sin(gamma) ) .* (sum((C7-IJ).^2,2)./sum((IJ-PX).^2,2)) , [1 3]) .* (IJ-PX);
        CJC = C7 + repmat( 0.55 .* ( sum((C7-IJ).^2,2) ./ sum((INTERSECT-C7).^2,2) ) , [1 3]) .* ( INTERSECT - C7 );
    elseif strcmpi(sex,'f')
        alpha = 14*pi/180;
        beta = acos( sum( (C7-IJ).*(IJ-PX),2 ) ./ ( sqrt( sum((C7-IJ).^2,2) ) .* sqrt( sum((IJ-PX).^2,2) ) ) );
        gamma = pi - alpha - beta;
        INTERSECT = IJ + repmat( ( sin(alpha)./sin(gamma) ) .* (sum((C7-IJ).^2,2)./sum((IJ-PX).^2,2)) , [1 3]) .* (IJ-PX);
        CJC = C7 + repmat( 0.53 .* ( sum((C7-IJ).^2,2) ./ sum((INTERSECT-C7).^2,2) ) , [1 3]) .* ( INTERSECT - C7 );
    end
    
    % y-axis : from LJC to CJC
    yvec = CJC - LJC;
    yvecNorm = yvec ./ repmat( sqrt(sum(yvec.^2,2)) , [1 3] );
    
    % z-axis : normal to plane LJC, CJC and IJ, pointing to the right
    zvec = cross(IJ - LJC,CJC - LJC);
    zvecNorm = zvec ./ repmat( sqrt(sum(zvec.^2,2)) , [1 3] );
    
    % x-axis : cross product y and z axes
    xvecNorm = cross(yvecNorm,zvecNorm);
    
    % Get rotation matrix to get from a segment coordinates to global coordinates
    seg2globR(:,n,1:3,1) = xvecNorm;
    seg2globR(:,n,1:3,2) = yvecNorm;
    seg2globR(:,n,1:3,3) = zvecNorm;
    
    % Origin (CJC)
    seg2globD(:,n,:) = CJC;
    
    % Joints
    m = m + 1;
    jointData(:,m,:) = CJC;
    jointDataLabel{m} = 'CJC';
    
    % Segment length
    segLen = sqrt( sum((CJC - LJC).^2,2) );
    
    % Inertia and COM local, mass
    if strcmpi(sex,'m')
        comLocal = [-0.036 -0.420 0]';
        iMat = [0.27 0.18 0.02 ; 0.18 0.25 -0.04 ; 0.02 -0.04 0.28];
        segMass(n) = 0.333 .* mass;
    elseif strcmpi(sex,'f')
        comLocal = [-0.016 -0.436 0]';
        iMat = [0.29 0.22 0.05 ; 0.22 0.27 -0.05 ; 0.05 -0.05 0.29];
        segMass(n) = 0.304 .* mass;
    end
    
    % Inertia and COM global
    for i = 1:size(seg2globR,1)
        comData(i,n,:) = reshape(seg2globD(i,n,:),[3 1]) + reshape(seg2globR(i,n,:,:),[3 3]) * (segLen(i).*comLocal);
        
        % For alternative origin of the trunk inertia matrix: LJC
        % (disable if unused)
%         LJClocal = [0 -sqrt(sum(yvec(i,:).^2)) 0].';
%         parAxShft = eye(3).*norm(LJClocal).^2 - (LJClocal)*(LJClocal.');
%         segInertiaO(i,n,:,:) = reshape(seg2globR(i,n,:,:),[3 3]) * ( segMass(n).* ((iMat*segLen(i)).^2 + parAxShft) ) * reshape(seg2globR(i,n,:,:),[3 3])'; % About COM        
        
        segInertiaO(i,n,:,:) = reshape(seg2globR(i,n,:,:),[3 3]) * ( segMass(n).* (iMat*segLen(i)).^2 ) * reshape(seg2globR(i,n,:,:),[3 3])'; % About origin
        parAxShft = eye(3).*norm(segLen(i).*comLocal).^2 - (segLen(i).*comLocal)*(segLen(i).*comLocal.');
        segInertia(i,n,:,:) = reshape(seg2globR(i,n,:,:),[3 3]) * ( segMass(n).* ((iMat*segLen(i)).^2 + parAxShft) ) * reshape(seg2globR(i,n,:,:),[3 3])'; % About COM

    end
    
    % Obtain the names for COMLabelCell, using the cluster names
    comDataLabel{n} = labelallC{strcmpi('C7',labelall)};

end

%% HEAD

if isHEAD
    n = n + 1;
    
    % NOTE : requires the CJC from the TRUNK as origin
    
    % y-axis : from CJC to HV
    yvec = HV - CJC;
    yvecNorm = yvec ./ repmat( sqrt(sum(yvec.^2,2)) , [1 3] );
    
    % z-axis : normal to plane HV, CJC and SEL, pointing lateral (right)
    zvec = cross(HV - CJC , SEL - CJC);
    zvecNorm = zvec ./ repmat( sqrt(sum(zvec.^2,2)) , [1 3] );
    
    % x-axis : cross product
    xvecNorm = cross(yvecNorm,zvecNorm);
    
    % Get rotation matrix to get from a segment coordinates to global coordinates
    seg2globR(:,n,1:3,1) = xvecNorm;
    seg2globR(:,n,1:3,2) = yvecNorm;
    seg2globR(:,n,1:3,3) = zvecNorm;
    
    % Origin (CJC)
    seg2globD(:,n,:) = CJC;
    
    % Segment length
    segLen = sqrt( sum((CJC - HV).^2,2) );
    
    % Inertia and COM local
    if strcmpi(sex,'m')
        comLocal = [-0.062 0.555 0]';
        iMat = [0.31 -0.09 -0.02 ; -0.09 0.25 0.03 ; -0.02 0.03 0.33];
        segMass(n) = 0.067 .* mass;
    elseif strcmpi(sex,'f')
        comLocal = [-0.070 0.597 0]';
        iMat = [0.32 -0.06 0.01 ; -0.06 0.27 -0.01 ; 0.01 -0.01 0.34];
        segMass(n) = 0.067 .* mass;
    end
    
    % Inertia and COM global
    for i = 1:size(seg2globR,1)
        comData(i,n,:) = reshape(seg2globD(i,n,:),[3 1]) + reshape(seg2globR(i,n,:,:),[3 3]) * (segLen(i).*comLocal);
        
        segInertiaO(i,n,:,:) = reshape(seg2globR(i,n,:,:),[3 3]) * ( segMass(n).* (iMat*segLen(i)).^2 ) * reshape(seg2globR(i,n,:,:),[3 3])'; % About origin
        parAxShft = eye(3).*norm(segLen(i).*comLocal).^2 - (segLen(i).*comLocal)*(segLen(i).*comLocal.');
        segInertia(i,n,:,:) = reshape(seg2globR(i,n,:,:),[3 3]) * ( segMass(n).* ((iMat*segLen(i)).^2 + parAxShft) ) * reshape(seg2globR(i,n,:,:),[3 3])'; % About COM
    end
    
    % Obtain the names for COMLabelCell, using the cluster names
    comDataLabel{n} = labelallC{strcmpi('HV',labelall)};

end


%% LEFT SIDE
%% LFOOT
 
if isLFOOT
    n = n + 1;
    
    % Get AJC
    AJCL = (MLL + MML)./2;
    
    % x-axis : from CAL to midpoint CM1L and CM5L
    xvec = (CM1L+CM5L)./2 - CAL; 
    xvecNorm = xvec ./ repmat( sqrt(sum(xvec.^2,2)) , [1 3]) ;
    
    % y-axis : normal to the plane containing CAL, CM1L and CM5L, pointing cranial
    yvec = cross(CM1L - CAL , CM5L - CAL);
    yvecNorm = yvec ./ repmat( sqrt(sum(yvec.^2,2)) , [1 3]);
    
    % z-axis : cross product of x and y axes
    zvecNorm = cross(xvecNorm,yvecNorm);
    
    % Get rotation matrix to get from a segment coordinates to global coordinates
    seg2globR(:,n,1:3,1) = xvecNorm;
    seg2globR(:,n,1:3,2) = yvecNorm;
    seg2globR(:,n,1:3,3) = zvecNorm;
    
    % Origin (AJC)
    seg2globD(:,n,:) = AJCL;
    
    % Joints
    m = m + 1;
    jointData(:,m,:) = AJCL;
    jointDataLabel{m} = 'AJCL';
    
    % Segment length
    segLen = sqrt( sum((AJCL - (CM1L+CM5L)./2).^2,2) );
    
    % Inertia and COM local
    if strcmpi(sex,'m')
        comLocal = [0.382 -0.151 0.026]';
        iMat = [0.17 0.13 -0.08 ; 0.13 0.27 0 ; -0.08 0 0.36];
        segMass(n) = 0.012 .* mass;
    elseif strcmpi(sex,'f')
        comLocal = [0.270 -0.218 0.039]';
        iMat = [0.17 -0.10 0.06 ; -0.10 0.36 -0.04 ; 0.06 -0.04 0.35];
        segMass(n) = 0.010 .* mass;
    end

    % Inertia and COM global
    for i = 1:size(seg2globR,1)
        comData(i,n,:) = reshape(seg2globD(i,n,:),[3 1]) + reshape(seg2globR(i,n,:,:),[3 3]) * (segLen(i).*comLocal);
        
        segInertiaO(i,n,:,:) = reshape(seg2globR(i,n,:,:),[3 3]) * ( segMass(n).* (iMat*segLen(i)).^2 ) * reshape(seg2globR(i,n,:,:),[3 3])'; % About origin
        parAxShft = eye(3).*norm(segLen(i).*comLocal).^2 - (segLen(i).*comLocal)*(segLen(i).*comLocal.');
        segInertia(i,n,:,:) = reshape(seg2globR(i,n,:,:),[3 3]) * ( segMass(n).* ((iMat*segLen(i)).^2 + parAxShft) ) * reshape(seg2globR(i,n,:,:),[3 3])'; % About COM
    end
    
    % Obtain the names for COMLabelCell, using the cluster names
    comDataLabel{n} = labelallC{strcmpi('CAL',labelall)};
    
end


%% LLOWER

if isLLOWER
    n = n + 1;
    
    % Find AJC
    AJCL = (MLL + MML)./2;
    
    % Find KJC (origin)
    KJCL = (ELL + EML)./2;
    
    % y-axis : from AJC to KJC
    yvec = KJCL - AJCL;
    yvecNorm = yvec ./ repmat( sqrt(sum(yvec.^2,2)) , [1 3]);
    
    % x-axis : normal to the plane AJC-KJC-CFL, pointing anterior
    xvec = cross(AJCL - KJCL , CFL - KJCL);
    xvecNorm = xvec ./ repmat( sqrt(sum(xvec.^2,2)) , [1 3]);
    
    % z-axis : cross product of x and y axes
    zvecNorm = cross(xvecNorm,yvecNorm);
    
    % Get rotation matrix to get from a segment coordinates to global coordinates
    seg2globR(:,n,1:3,1) = xvecNorm;
    seg2globR(:,n,1:3,2) = yvecNorm;
    seg2globR(:,n,1:3,3) = zvecNorm;
    
    % Origin
    seg2globD(:,n,:) = KJCL;
    
    % Joints
    if ~any(strcmpi('AJCL',jointDataLabel)) % Check if not already assigned in LFOOT
        m = m + 1;
        jointData(:,m,:) = AJCL;
        jointDataLabel{m} = 'AJCL';
    end
    
    m = m + 1;
    jointData(:,m,:) = KJCL;
    jointDataLabel{m} = 'KJCL';
    
    % Segment length
    segLen = sqrt( sum((KJCL - AJCL).^2,2) );

    % Inertia and COM local
    if strcmpi(sex,'m')
        comLocal = [-0.048 -0.410 0]';
        iMat = [0.28 -0.04 -0.02 ; -0.04 0.10 0.05 ; -0.02 0.05 0.28];
        segMass(n) = 0.048 .* mass;
    elseif strcmpi(sex,'f')
        comLocal = [-0.049 -0.404 0.031]';
        iMat = [0.28 0.02 0.01 ; 0.02 0.10 0.06 ; 0.01 0.06 0.28];
        segMass(n) = 0.045 .* mass;
    end
    
    % Inertia and COM global
    for i = 1:size(seg2globR,1)
        comData(i,n,:) = reshape(seg2globD(i,n,:),[3 1]) + reshape(seg2globR(i,n,:,:),[3 3]) * (segLen(i).*comLocal);
        
        segInertiaO(i,n,:,:) = reshape(seg2globR(i,n,:,:),[3 3]) * ( segMass(n).* (iMat*segLen(i)).^2 ) * reshape(seg2globR(i,n,:,:),[3 3])'; % About origin
        parAxShft = eye(3).*norm(segLen(i).*comLocal).^2 - (segLen(i).*comLocal)*(segLen(i).*comLocal.');
        segInertia(i,n,:,:) = reshape(seg2globR(i,n,:,:),[3 3]) * ( segMass(n).* ((iMat*segLen(i)).^2 + parAxShft) ) * reshape(seg2globR(i,n,:,:),[3 3])'; % About COM
    end
    
    % Obtain the names for COMLabelCell, using the cluster names
    comDataLabel{n} = labelallC{strcmpi('CFL',labelall)};
    
end


%% LUPPER

if isLUPPER
    n = n + 1;
    % NOTE : the hip joint centers must be found in the PELVIS section
    
    % Find KJC (origin)
    KJCL = (ELL + EML)./2;
    
    % y-axis : from KJC to HJC
    yvec = HJCL - KJCL;
    yvecNorm = yvec ./ repmat( sqrt(sum(yvec.^2,2)) , [1 3]);
    
    % x-axis : normal to plane containing HJCL, ELL, EML, pointing anterior
    xvec = cross( EML - HJCL , ELL - HJCL );
    xvecNorm = xvec ./ repmat( sqrt(sum(xvec.^2,2)) , [1 3]);
    
    % z-axis : cross product of x and y axes
    zvecNorm = cross(xvecNorm,yvecNorm);
    
    % Get rotation matrix to get from a segment coordinates to global coordinates
    seg2globR(:,n,1:3,1) = xvecNorm;
    seg2globR(:,n,1:3,2) = yvecNorm;
    seg2globR(:,n,1:3,3) = zvecNorm;
    
    % Origin
    seg2globD(:,n,:) = HJCL;
    
    % Joints
    if ~any(strcmpi('KJCL',jointDataLabel))
        m = m + 1;
        jointData(:,m,:) = KJCL;
        jointDataLabel{m} = 'KJCL'; 
    end
    
    % Segment length
    segLen = sqrt( sum((KJCL - HJCL).^2,2) );
    
    % Inertia and COM local
    if strcmpi(sex,'m')
        comLocal = [-0.041 -0.429 0.033]';
        iMat = [0.29 0.07 -0.02 ; 0.07 0.15 -0.07 ; -0.02 -0.07 0.30];
        segMass(n) = 0.123 .* mass;
    elseif strcmpi(sex,'f')
        comLocal = [0.077 -0.377 0]';
        iMat = [0.31 0.07 -0.02 ; 0.07 0.19 -0.07 ; -0.02 -0.07 0.32];
        segMass(n) = 0.146 .* mass;
    end
    
    % Inertia and COM global
    for i = 1:size(seg2globR,1)
        comData(i,n,:) = reshape(seg2globD(i,n,:),[3 1]) + reshape(seg2globR(i,n,:,:),[3 3]) * (segLen(i).*comLocal);
        
        segInertiaO(i,n,:,:) = reshape(seg2globR(i,n,:,:),[3 3]) * ( segMass(n).* (iMat*segLen(i)).^2 ) * reshape(seg2globR(i,n,:,:),[3 3])'; % About origin
        parAxShft = eye(3).*norm(segLen(i).*comLocal).^2 - (segLen(i).*comLocal)*(segLen(i).*comLocal.');
        segInertia(i,n,:,:) = reshape(seg2globR(i,n,:,:),[3 3]) * ( segMass(n).* ((iMat*segLen(i)).^2 + parAxShft) ) * reshape(seg2globR(i,n,:,:),[3 3])'; % About COM
    end
    
    % Obtain the names for COMLabelCell, using the cluster names
    comDataLabel{n} = labelallC{strcmpi('ELL',labelall)};
    
end


%% RIGHT SIDE
%% RFOOT
 
if isRFOOT
    n = n + 1;
    
    % Get AJC
    AJCR = (MLR + MMR)./2;
    
    % x-axis : from CAL to midpoint CM1L and CM5L
    xvec = (CM1R+CM5R)./2 - CAR; 
    xvecNorm = xvec ./ repmat( sqrt(sum(xvec.^2,2)) , [1 3]) ;
    
    % y-axis : normal to the plane containing CAL, CM1L and CM5L, pointing cranial
    yvec = cross(CM5R - CAR , CM1R - CAR);
    yvecNorm = yvec ./ repmat( sqrt(sum(yvec.^2,2)) , [1 3]);
    
    % z-axis : cross product of x and y axes
    zvecNorm = cross(xvecNorm,yvecNorm);
    
    % Get rotation matrix to get from a segment coordinates to global coordinates
    seg2globR(:,n,1:3,1) = xvecNorm;
    seg2globR(:,n,1:3,2) = yvecNorm;
    seg2globR(:,n,1:3,3) = zvecNorm;
    
    % Origin (AJC)
    seg2globD(:,n,:) = AJCR;
    
    % Joints
    m = m + 1;
    jointData(:,m,:) = AJCR;
    jointDataLabel{m} = 'AJCR';
    
    % Segment length
    segLen = sqrt( sum((AJCR - (CM1R+CM5R)./2).^2,2) );
    
    % Inertia and COM local
    if strcmpi(sex,'m')
        comLocal = [0.382 -0.151 0.026]';
        iMat = [0.17 0.13 -0.08 ; 0.13 0.27 0 ; -0.08 0 0.36];
        segMass(n) = 0.012 .* mass;
    elseif strcmpi(sex,'f')
        comLocal = [0.270 -0.218 0.039]';
        iMat = [0.17 -0.10 0.06 ; -0.10 0.36 -0.04 ; 0.06 -0.04 0.35];
        segMass(n) = 0.010 .* mass;
    end
    
    % Inertia and COM global
    for i = 1:size(seg2globR,1)
        comData(i,n,:) = reshape(seg2globD(i,n,:),[3 1]) + reshape(seg2globR(i,n,:,:),[3 3]) * (segLen(i).*comLocal);
        
        segInertiaO(i,n,:,:) = reshape(seg2globR(i,n,:,:),[3 3]) * ( segMass(n).* (iMat*segLen(i)).^2 ) * reshape(seg2globR(i,n,:,:),[3 3])'; % About origin
        parAxShft = eye(3).*norm(segLen(i).*comLocal).^2 - (segLen(i).*comLocal)*(segLen(i).*comLocal.');
        segInertia(i,n,:,:) = reshape(seg2globR(i,n,:,:),[3 3]) * ( segMass(n).* ((iMat*segLen(i)).^2 + parAxShft) ) * reshape(seg2globR(i,n,:,:),[3 3])'; % About COM
    end
    
    % Obtain the names for COMLabelCell, using the cluster names
    comDataLabel{n} = labelallC{strcmpi('CAR',labelall)};
    
end

    
%% RLOWER

if isRLOWER
    n = n + 1;
    
    % Find AJC
    AJCR = (MLR + MMR)./2;
    
    % Find KJC (origin)
    KJCR = (ELR + EMR)./2;
    
    % y-axis : from AJC to KJC
    yvec = KJCR - AJCR;
    yvecNorm = yvec ./ repmat( sqrt(sum(yvec.^2,2)) , [1 3]);
    
    % x-axis : normal to the plane AJC-KJC-CFL, pointing anterior
    xvec = cross(CFR - KJCR , AJCR - KJCR);
    xvecNorm = xvec ./ repmat( sqrt(sum(xvec.^2,2)) , [1 3]);
    
    % z-axis : cross product of x and y axes
    zvecNorm = cross(xvecNorm,yvecNorm);
    
    % Get rotation matrix to get from a segment coordinates to global coordinates
    seg2globR(:,n,1:3,1) = xvecNorm;
    seg2globR(:,n,1:3,2) = yvecNorm;
    seg2globR(:,n,1:3,3) = zvecNorm;
    
    % Origin
    seg2globD(:,n,:) = KJCR;
    
    % Joints
    if ~any(strcmpi('AJCR',jointDataLabel)) % Check if not already assigned in LFOOT
        m = m + 1;
        jointData(:,m,:) = AJCR;
        jointDataLabel{m} = 'AJCR';
    end
    
    m = m + 1;
    jointData(:,m,:) = KJCR;
    jointDataLabel{m} = 'KJCR';
    
    % Segment length
    segLen = sqrt( sum((KJCR - AJCR).^2,2) );
    
    % Inertia and COM local
    if strcmpi(sex,'m')
        comLocal = [-0.048 -0.410 0]';
        iMat = [0.28 -0.04 -0.02 ; -0.04 0.10 0.05 ; -0.02 0.05 0.28];
        segMass(n) = 0.048 .* mass;
    elseif strcmpi(sex,'f')
        comLocal = [-0.049 -0.404 0.031]';
        iMat = [0.28 0.02 0.01 ; 0.02 0.10 0.06 ; 0.01 0.06 0.28];
        segMass(n) = 0.045 .* mass;
    end
    
    % Inertia and COM global
    for i = 1:size(seg2globR,1)
        comData(i,n,:) = reshape(seg2globD(i,n,:),[3 1]) + reshape(seg2globR(i,n,:,:),[3 3]) * (segLen(i).*comLocal);
        
        segInertiaO(i,n,:,:) = reshape(seg2globR(i,n,:,:),[3 3]) * ( segMass(n).* (iMat*segLen(i)).^2 ) * reshape(seg2globR(i,n,:,:),[3 3])'; % About origin
        parAxShft = eye(3).*norm(segLen(i).*comLocal).^2 - (segLen(i).*comLocal)*(segLen(i).*comLocal.');
        segInertia(i,n,:,:) = reshape(seg2globR(i,n,:,:),[3 3]) * ( segMass(n).* ((iMat*segLen(i)).^2 + parAxShft) ) * reshape(seg2globR(i,n,:,:),[3 3])'; % About COM
    end
    
    % Obtain the names for COMLabelCell, using the cluster names
    comDataLabel{n} = labelallC{strcmpi('CFR',labelall)};
    
end

%% RUPPER

if isRUPPER
    n = n + 1;
    % NOTE : the hip joint centers must be found in the PELVIS section
    
    % Find KJC (origin)
    KJCR = (ELR + EMR)./2;
    
    % y-axis : from KJC to HJC
    yvec = HJCR - KJCR;
    yvecNorm = yvec ./ repmat( sqrt(sum(yvec.^2,2)) , [1 3]);
    
    % x-axis : normal to plane containing HJCR, ELR, EMR, pointing anterior
    xvec = cross( ELR - HJCR , EMR - HJCR );
    xvecNorm = xvec ./ repmat( sqrt(sum(xvec.^2,2)) , [1 3]);
    
    % z-axis : cross product of x and y axes
    zvecNorm = cross(xvecNorm,yvecNorm);
    
    % Get rotation matrix to get from a segment coordinates to global coordinates
    seg2globR(:,n,1:3,1) = xvecNorm;
    seg2globR(:,n,1:3,2) = yvecNorm;
    seg2globR(:,n,1:3,3) = zvecNorm;
    
    % Origin
    seg2globD(:,n,:) = HJCR;
    
    % Joints
    if ~any(strcmpi('KJCR',jointDataLabel))
        m = m + 1;
        jointData(:,m,:) = KJCR;
        jointDataLabel{m} = 'KJCR'; 
    end
    
    % Segment length
    segLen = sqrt( sum((KJCR - HJCR).^2,2) );
    
    % Inertia and COM local
    if strcmpi(sex,'m')
        comLocal = [-0.041 -0.429 0.033]';
        iMat = [0.29 0.07 -0.02 ; 0.07 0.15 -0.07 ; -0.02 -0.07 0.30];
        segMass(n) = 0.123 .* mass;
    elseif strcmpi(sex,'f')
        comLocal = [0.077 -0.377 0]';
        iMat = [0.31 0.07 -0.02 ; 0.07 0.19 -0.07 ; -0.02 -0.07 0.32];
        segMass(n) = 0.146 .* mass;
    end
    
    % Inertia and COM global
    for i = 1:size(seg2globR,1)
        comData(i,n,:) = reshape(seg2globD(i,n,:),[3 1]) + reshape(seg2globR(i,n,:,:),[3 3]) * (segLen(i).*comLocal);
        
        segInertiaO(i,n,:,:) = reshape(seg2globR(i,n,:,:),[3 3]) * ( segMass(n).* (iMat*segLen(i)).^2 ) * reshape(seg2globR(i,n,:,:),[3 3])'; % About origin
        parAxShft = eye(3).*norm(segLen(i).*comLocal).^2 - (segLen(i).*comLocal)*(segLen(i).*comLocal.');
        segInertia(i,n,:,:) = reshape(seg2globR(i,n,:,:),[3 3]) * ( segMass(n).* ((iMat*segLen(i)).^2 + parAxShft) ) * reshape(seg2globR(i,n,:,:),[3 3])'; % About COM
    end
    
    % Obtain the names for COMLabelCell, using the cluster names
    comDataLabel{n} = labelallC{strcmpi('ELR',labelall)};
    
end


%% Weighted COM of the available segments
if nseg > 0
    n = n + 1;

    comData(:,n,:) = sum( comData(:,1:nseg,:) .* repmat(segMass,[size(comData,1) 1 size(comData,3)]) ,2) ./ sum(segMass);
    segMass(n) = sum(segMass);
    comDataLabel{n} = 'comALL';

end

%% Assign data

C3Ddata.Marker.JointData = jointData;
C3Ddata.Marker.JointDataLabel = jointDataLabel;

C3Ddata.Marker.COMData = comData;       % Global COM positions
C3Ddata.Marker.COMDataLabel = comDataLabel;

C3Ddata.Marker.SegInertiaO = segInertiaO; % Global segment inertias about origin. The comDataLabel applies to these too
C3Ddata.Marker.SegInertia = segInertia; % Global segment inertias about COM. The comDataLabel applies to these too
C3Ddata.Marker.SegMass = segMass;

C3Ddata.Marker.Seg2GlobR = seg2globR;  % Rotation from segment coordinates to global coordinates
C3Ddata.Marker.Seg2GlobD = seg2globD;  % Displacement from segment origin (0,0,0) to that in global coordinates

C3Ddata.Marker = orderfields(C3Ddata.Marker);


end