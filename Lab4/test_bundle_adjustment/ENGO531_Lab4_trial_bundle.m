%% ENGO 531 Lab 4 -- TRIAL bundle adjustment (NO DATA SNOOPING)
% Standalone nonlinear adjustment; requires Optimization Toolbox (lsqnonlin).
% IMAGE 1 EOP is fixed to establish datum; F1--F6 provides scale.
% IOPs are fixed. All XYZ other than camera centres are in METRES internally.
% No target coordinates are fabricated: starting XYZ are triangulated from images.
% IMPORTANT: confirm camera-angle and collinearity conventions with your lab.
clear; clc; close all;

%% SETTINGS
phoFile = 'calibration(2).PHO'; % also supports calibration_points(2).csv
imageWidth = 5712; imageHeight = 4284;
minImages = 3;
allowTwoImageScaleEndpoint = true; % TRIAL ONLY: F1 has only TWO images
sdImagePx = 3.2;
sdScaleM = 0.005; % provisional tape measurement uncertainty; change if known
scaleFrom = 'F1'; scaleTo = 'F6'; scaleLengthM = 2.478;
runBothIOPs = false; % true to compare nominal and calibrated IOPs
calibratedIOP = [2834 2156.7 3995.9]; % xp yp c (pixels)
nominalIOP = [imageWidth/2 imageHeight/2 4109.10]; % verify c from Lab 3

% INPUT COLUMN ORDER from your field notes: X, Y, Z, phi, kappa, omega.
fieldNotes = [ ...
 1 -2213  1711 -645    90    0 90;
 2 -2279   561 -495    75    0 90;
 3 -1619 -1540 -495    45    0 90;
 4  -440 -2027 -495    25    0 90;
 5  1239 -2545 -495     0    0 90;
 6  3486 -1652 -495   -30    0 90;
 7  4368   -97 -495   -60    0 90;
 8  4901   411 -495   -90    0 90;
 9  4455  2686 -495  -120    0 90];

%% LOAD ALPHANUMERIC IDs AND VALIDATE OVERLAP
if ~isfile(phoFile)
    [fn,fd] = uigetfile({'*.PHO;*.pho;*.csv','Image observations'},'Choose PHO/CSV');
    if isequal(fn,0), return; end
    phoFile = fullfile(fd,fn);
end
[~,~,extn] = fileparts(phoFile);
if strcmpi(extn,'.csv')
    raw = readtable(phoFile,'TextType','string');
    assert(all(ismember({'point_id','image_id','x_pix','y_pix'},raw.Properties.VariableNames)), ...
        'CSV must have point_id,image_id,x_pix,y_pix');
    pointName = string(raw.point_id);
    imageID = double(raw.image_id);
    measuredXY = [double(raw.x_pix) double(raw.y_pix)];
else
    fid = fopen(phoFile,'r'); assert(fid>=0,'Cannot open observations');
    C = textscan(fid,'%s %f %f %f'); fclose(fid);
    pointName = string(C{1}); imageID = C{2}; measuredXY = [C{3} C{4}];
end
assert(all(isfinite(imageID)) && all(isfinite(measuredXY(:))), 'Non-finite observations');
assert(all(measuredXY(:,1)>=0 & measuredXY(:,1)<=imageWidth & ...
           measuredXY(:,2)>=0 & measuredXY(:,2)<=imageHeight), 'Image coordinates outside dimensions');
assert(numel(unique(pointName+"_"+string(imageID)))==numel(imageID), 'Duplicate point/image observation');

% Include only actually observed images -- image 6 has no measurements.
activeImages = intersect(fieldNotes(:,1),unique(imageID),'stable');
missingImages = setdiff(fieldNotes(:,1),activeImages);
fprintf('Images with observations: %s\n',mat2str(activeImages'));
if ~isempty(missingImages)
    fprintf('Images without observations (omitted): %s\n',mat2str(missingImages'));
end

[allNames,~,g] = unique(pointName,'stable');
counts = accumarray(g,1);
ct = table(allNames,counts,'VariableNames',{'PointID','ImageCount'});
disp('Initial overlap:'); disp(ct);
keepNames = allNames(counts>=minImages);
scaleNames = string({scaleFrom;scaleTo});
for k=1:numel(scaleNames)
    if ~ismember(scaleNames(k),keepNames)
        idx = find(allNames==scaleNames(k),1);
        if isempty(idx), error('Scale endpoint %s has no measurements',scaleNames(k)); end
        if counts(idx)==2 && allowTwoImageScaleEndpoint
            warning('TRIAL EXCEPTION: scale point %s occurs in only 2 images.',scaleNames(k));
            keepNames(end+1) = scaleNames(k); %#ok<SAGROW>
        else
            error(['Scale endpoint %s appears in %d images. Add measurements to reach %d, ' ...
                'or change scale endpoints to a separately measured pair.'], ...
                scaleNames(k),counts(idx),minImages);
        end
    end
end
keep = ismember(pointName,keepNames) & ismember(imageID,activeImages);
pointName = pointName(keep); imageID = imageID(keep); measuredXY = measuredXY(keep,:);
[pointIDs,~,pointIndex] = unique(pointName,'stable');
[tf,imageIndex] = ismember(imageID,activeImages);
assert(all(tf),'Internal image index failure');
nImages = numel(activeImages); nPoints = numel(pointIDs); nObs = numel(imageID);
fprintf('Using %d images, %d points, %d image-point measurements.\n',nImages,nPoints,nObs);
fprintf('All non-exception targets have at least %d images.\n',minImages);
assert(nImages>=2 && nPoints>=3,'Insufficient data');
[~,f1] = ismember(string(scaleFrom),pointIDs);
[~,f6] = ismember(string(scaleTo),pointIDs);
assert(f1>0 && f6>0,'Scale endpoints missing after filter');

% All images have EOP starting values. Convert mm to m; order omega,phi,kappa.
[tf,rows] = ismember(activeImages,fieldNotes(:,1)); assert(all(tf));
F = fieldNotes(rows,:);
EOP0 = [F(:,2:4)/1000 F(:,7) F(:,5) F(:,6)];
% Fix all six parameters of the first active image (Image 1).
% This is an exact constraint for a convergence trial, not a weighted observation.
fixedEOP = EOP0(1,:);

if runBothIOPs
    IOPsets = [calibratedIOP; nominalIOP];
    names = {'calibrated','nominal'};
else
    IOPsets = calibratedIOP;
    names = {'calibrated'};
end

for which = 1:size(IOPsets,1)
    IOP = IOPsets(which,:);
    fprintf('\n========== %s IOP: xp %.3f, yp %.3f, c %.3f ==========\n', ...
        names{which},IOP);
    % Ray-intersection initial point positions. All 3D point coordinates
    % are estimated from the supplied image measurements and approximate EOPs.
    XYZ0 = zeros(nPoints,3);
    for j = 1:nPoints
        obsIdx = find(pointIndex==j);
        A = zeros(3); b = zeros(3,1);
        for h = obsIdx'
            ii = imageIndex(h);
            C = EOP0(ii,1:3)';
            R = rotationWPK(deg2rad(EOP0(ii,4:6)));
            rayCamera = [measuredXY(h,1)-IOP(1); measuredXY(h,2)-IOP(2); -IOP(3)];
            rayWorld = R' * rayCamera;
            rayWorld = rayWorld/norm(rayWorld);
            M = eye(3)-rayWorld*rayWorld';
            A = A+M; b = b+M*C;
        end
        if rcond(A)<1e-12, error('Poor triangulation geometry for point %s',pointIDs(j)); end
        XYZ0(j,:) = (A\b)';
    end
    fprintf('Initial F1--F6 triangulated distance: %.4f m; measured %.4f m\n', ...
        norm(XYZ0(f1,:)-XYZ0(f6,:)),scaleLengthM);

    % Optimize EOPs for images 2..N, and all object point XYZ.
    % Positions in metres, angles in radians, XYZ in metres.
    eopParams = [EOP0(:,1:3) deg2rad(EOP0(:,4:6))];
    fixedParams = eopParams(1,:);
    x0 = [reshape(eopParams(2:end,:)',[],1); reshape(XYZ0',[],1)];
    residual = @(x) bundleResidual(x,fixedParams,nImages,nPoints, ...
                 pointIndex,imageIndex,measuredXY,IOP, ...
                 f1,f6,scaleLengthM,sdImagePx,sdScaleM);
    if exist('lsqnonlin','file')~=2
        error(['This standalone trial script needs Optimization Toolbox (lsqnonlin). ' ...
            'The course full_matlab_bundle.p is another option if you have it.']);
    end
    opts = optimoptions('lsqnonlin','Display','iter', ...
        'MaxIterations',150,'MaxFunctionEvaluations',30000, ...
        'FunctionTolerance',1e-10,'StepTolerance',1e-10, ...
        'FiniteDifferenceType','central');
    [xhat,resnorm,weightedResidual,exitflag,output] = lsqnonlin(residual,x0,[],[],opts);
    [eopFinal,XYZfinal] = unpack(xhat,fixedParams,nImages,nPoints);
    imageResidual = weightedResidual(1:2*nObs)*sdImagePx;
    rmsImagePx = sqrt(mean(imageResidual.^2));
    finalDistance = norm(XYZfinal(f1,:)-XYZfinal(f6,:));
    dof = 2*nObs+1-numel(xhat);
    fprintf('\nRESULT: exitflag %d; iterations %d; weighted SSE %.4f\n', ...
        exitflag,output.iterations,resnorm);
    fprintf('RMS %.4f px; F1--F6 %.4f m (tape %.4f m); DOF %d\n', ...
        rmsImagePx,finalDistance,scaleLengthM,dof);
    if exitflag<=0, warning('Solver did NOT report convergence. Inspect geometry and initialization.'); end
    if dof<=0, warning('Degrees of freedom <= 0; inspect overlap/constraints.'); end
    if rmsImagePx>10, warning('High residual RMS: verify EOP angle and image coordinate conventions!'); end
    outPoints = table(pointIDs,XYZfinal(:,1)*1000,XYZfinal(:,2)*1000,XYZfinal(:,3)*1000, ...
        'VariableNames',{'PointID','X_mm','Y_mm','Z_mm'});
    outEOP = table(activeImages,eopFinal(:,1)*1000,eopFinal(:,2)*1000,eopFinal(:,3)*1000, ...
         rad2deg(eopFinal(:,4)),rad2deg(eopFinal(:,5)),rad2deg(eopFinal(:,6)), ...
         'VariableNames',{'ImageID','X_mm','Y_mm','Z_mm','Omega_deg','Phi_deg','Kappa_deg'});
    disp(outEOP); disp(outPoints);
    folder = fileparts(whichFileAbsolute(phoFile));
    writetable(outPoints,fullfile(folder,['Lab4_trial_points_' names{which} '.csv']));
    writetable(outEOP,fullfile(folder,['Lab4_trial_EOP_' names{which} '.csv']));
    save(fullfile(folder,['Lab4_trial_results_' names{which} '.mat']), ...
        'outPoints','outEOP','imageResidual','exitflag','output', ...
        'rmsImagePx','finalDistance','dof','pointName','imageID','IOP');
end

%% LOCAL FUNCTIONS
function R = rotationWPK(wpk)
 w=wpk(1); p=wpk(2); k=wpk(3);
 Rx=[1 0 0;0 cos(w) -sin(w);0 sin(w) cos(w)];
 Ry=[cos(p) 0 sin(p);0 1 0;-sin(p) 0 cos(p)];
 Rz=[cos(k) -sin(k) 0;sin(k) cos(k) 0;0 0 1];
 R=Rz*Ry*Rx;
end
function [E,P] = unpack(x,fixed,nImages,nPoints)
 n=6*(nImages-1);
 E=[fixed;reshape(x(1:n),6,[])'];
 P=reshape(x(n+1:n+3*nPoints),3,[])';
end
function r = bundleResidual(x,fixed,nImages,nPoints,pointIndex,imageIndex,measuredXY,IOP,p1,p2,d,sdImage,sdDist)
 [E,P]=unpack(x,fixed,nImages,nPoints);
 n=numel(pointIndex); predicted=zeros(n,2);
 for i=1:nImages
     ii=find(imageIndex==i);
     if isempty(ii), continue; end
     R=rotationWPK(E(i,4:6));
     q=R*(P(pointIndex(ii),:)-E(i,1:3))';
     dz=q(3,:)';
     dz(abs(dz)<1e-7)=sign(dz(abs(dz)<1e-7)+eps)*1e-7;
     predicted(ii,1)=IOP(1)-IOP(3)*q(1,:)'./dz;
     predicted(ii,2)=IOP(2)-IOP(3)*q(2,:)'./dz;
 end
 r=[reshape(((measuredXY-predicted)/sdImage)',[],1);
    (norm(P(p1,:)-P(p2,:))-d)/sdDist];
end
function p = whichFileAbsolute(f)
 if isfile(f), p=which(f); else, p=''; end
 if isempty(p), p=fullfile(pwd,f); end
end
