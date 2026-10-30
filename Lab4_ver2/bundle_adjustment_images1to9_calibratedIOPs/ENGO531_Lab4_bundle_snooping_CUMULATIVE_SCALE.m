%% ENGO 531 Lab 4 - interactive data snooping + cumulative variance scaling
% Calls full_matlab_bundle.p; allows removal of target/image observations or whole targets.
% Uses Lab 1 mode 1: FIRST-image weighted EOP observations + two measured distances.
% The photogrammetric IOPs are fixed to Lab 3 estimates.
% Tie-point approximate XYZ are initialized from the photographed pod faces.
clear; clc; close all;

%% SETTINGS
imageSize = [5712 4284];      % width,height in pixels
IOP = [2834.0 2156.7 3995.9]; % xp,yp,c (pixels); fixed
% Change this between runs. Examples: [1 2 3 4 5 6], [1:7], [1:9].
% Datum: the FIRST selected image gets six weighted EOP observations.
% Include Image 1 in comparisons if you want to preserve the same datum.
selectedImages = [1:9];
minImagesPerPoint = 3;       % Set to 2 for additional troubleshooting
% Change to true only to allow F1/F6 with two images when minimum is 3.
allowTwoImageScaleEndpoints = false;
scaleFrom = "F1"; scaleTo = "F6";
measuredDistance_mm = 2471;
secondFrom = "R4"; secondTo = "R1";
secondDistance_mm = 2701;
sd_xy = 1.0;                  % image observations (pixels)
sd_dist = 5;                 % trial tape-measure SD (mm)
sd_firstEOP_XYZ = 1;         % first camera EOP weighted position SD (mm)
sd_firstEOP_wpk = 0.01;      % first camera angles SD (degrees)

% Original field-note EOPs with reversed phi: [image Xc Yc Zc phi kappa omega]
notes = [
 1 -2213.00  1711.00 -645.00  -90.000   0.000  90.000;
 2 -1946.93  -307.61 -697.22  -55.209   0.956  86.608;
 3  -845.00 -1543.16 -592.78  -34.096   4.934  86.519;
 4   378.74 -2012.41 -388.77  -20.438   6.661  83.917;
 5  1772.39 -1795.98 -202.55    4.932  10.906  81.931;
 6  3941.82  -933.02  170.01   44.496  26.373  67.234;
 7  4888.63   832.82  397.52   63.879  47.547  45.247;
 8  4794.84  2075.29  418.94   71.966  99.436  -5.958;
 9  4558.46  3371.02  421.54   55.442 150.123 -55.516];

%% INPUT FILES
[fn,fd] = uigetfile({'*.PHO;*.pho','Image point measurements (*.PHO)'}, ...
    'Select the UPDATED calibration(3).PHO with Image 6 observations');
if isequal(fn,0), return; end
phoFile = fullfile(fd,fn);
outDir = uigetdir(fd,'Select output folder for selected-image convergence trial');
if isequal(outDir,0), return; end
imageTag = strjoin(string(selectedImages),'_');
prefix = fullfile(outDir, ['engo531_lab4_images_' char(imageTag)]);

fid=fopen(phoFile,'r'); assert(fid>=0,'Cannot open PHO');
c=textscan(fid,'%s %f %f %f'); fclose(fid);
pointName=string(c{1}); imageID=c{2}; xy=[c{3},c{4}];
assert(~isempty(imageID),'Empty PHO file');
assert(all(isfinite([imageID;xy(:)])),'PHO contains missing/non-numeric data');
assert(all(xy(:,1)>=0 & xy(:,1)<=imageSize(1) & ...
           xy(:,2)>=0 & xy(:,2)<=imageSize(2)),'Image measurements outside format');
assert(numel(unique(pointName+"_"+string(imageID)))==numel(imageID), ...
    'Duplicate target/image measurement');
% Select images BEFORE counting target overlap.
selectedImages = selectedImages(:)';
assert(numel(selectedImages)>=2 && ...
    all(ismember(selectedImages,notes(:,1))) && ...
    numel(unique(selectedImages))==numel(selectedImages), ...
    'Choose at least two distinct image IDs between 1 and 9');
assert(isequal(selectedImages,sort(selectedImages)), ...
    'List selected images in ascending order');
assert(all(ismember(selectedImages,unique(imageID))), ...
    'Selected images are missing from PHO');
activeImages=selectedImages(:);
keepImage = ismember(imageID,activeImages);
pointName=pointName(keepImage);
imageID=imageID(keepImage);
xy=xy(keepImage,:);
fprintf('Selected images: %s\n',mat2str(activeImages'));

%% RETAIN POINTS SEEN IN THE SELECTED SUBSET
[names0,~,g]=unique(pointName,'stable'); counts=accumarray(g,1);
disp(table(names0,counts,'VariableNames',{'Target','ImageCount'}));
keepNames=names0(counts>=minImagesPerPoint);
for endpoint=[scaleFrom;scaleTo]'
  j=find(names0==endpoint,1);
  if ~isempty(j) && counts(j)>=2 && allowTwoImageScaleEndpoints && ...
          counts(j)<minImagesPerPoint
      keepNames(end+1)=endpoint; %#ok<AGROW>
      warning('Temporary scale-point exception: %s has only %d images',endpoint,counts(j));
  end
end
sel=ismember(pointName,keepNames);
pointName=pointName(sel); imageID=imageID(sel); xy=xy(sel,:);
[names,~,pointIdx]=unique(pointName,'stable');
numPoints=numel(names); numImages=numel(activeImages); numPhoto=size(xy,1);
% Check connectivity after overlap filtering. This is necessary but does not
% guarantee good intersection geometry or a full-rank normal matrix.
visited=false(numImages,1); visited(1)=true;
changed=true;
while changed
    changed=false;
    for a=1:numImages
        if ~visited(a),continue;end
        observed=unique(pointName(imageID==activeImages(a)));
        for b=1:numImages
            if ~visited(b) && any(ismember(observed,pointName(imageID==activeImages(b))))
                visited(b)=true;changed=true;
            end
        end
    end
end
assert(all(visited), ...
    'Selected network is DISCONNECTED after filtering. Add bridging images or reduce minimum overlap.');

assert(all(ismember([scaleFrom;scaleTo],names)), ...
    'F1/F6 are missing from this image subset: cannot use the existing scale distance. Select more images.');
assert(all(ismember(activeImages,unique(imageID))), ...
    'Overlap filter removed every point from at least one selected image');

% Number every target uniquely for course P-code; keep reversible CSV mapping.
numericID=(1001:(1000+numPoints))';
idMap=table(names,numericID,'VariableNames',{'OriginalPointID','BundlePointID'});
writetable(idMap,[prefix '_target_ID_map.csv']);

%% INITIAL OBJECT XYZ FROM PHOTOGRAPHED FACE GEOMETRY (mm)
% Coordinate convention of prior successful trial:
%   Front face Y = 0; width 2478 mm, height 2213 mm (Z downward negative).
%   Left wall X approx -65, Y from 0 (front) to +2704 (back).
%   Right wall X approx +2543, Y from 0 (front) to +2704 (back).
% These are APPROXIMATE starting values, not control coordinates.
% Four visible corners per face establish projective image-to-face mappings.
leftH = faceHomography(["L1" "L4" "L5" "L8"], 1, ...
    pointNameFull(phoFile), imageIDsFull(phoFile), coordinatesFull(phoFile), ...
    [0 0;2704 0;2704 -2213;0 -2213]);
frontH = faceHomography(["F1" "F6" "F7" "F12"], 5, ...
    pointNameFull(phoFile), imageIDsFull(phoFile), coordinatesFull(phoFile), ...
    [0 0;2478 0;2478 -2213;0 -2213]);
rightH = faceHomography(["R4" "R1" "R8" "R5"], 8, ...
    pointNameFull(phoFile), imageIDsFull(phoFile), coordinatesFull(phoFile), ...
    [0 0;2704 0;2704 -2213;0 -2213]);

allNames=pointNameFull(phoFile);
allImages=imageIDsFull(phoFile);
allXY=coordinatesFull(phoFile);
XYZ0=zeros(numPoints,3);
for j=1:numPoints
    nm=names(j);
    if startsWith(nm,"L")
        im=1; H=leftH;
        % L12 is not visible in image 1; use a simple geometric estimate.
        if nm=="L12"
            XYZ0(j,:)=[-65 1352 -2213];
            continue;
        end
        q=find(allNames==nm & allImages==im,1);
        if isempty(q)
            im=3; q=find(allNames==nm & allImages==im,1);
            % Image 3 point on left face needs a homography for Image 3.
            H=faceHomography(["L1" "L4" "L5" "L8"],3,allNames,allImages,allXY, ...
               [0 0;2704 0;2704 -2213;0 -2213]);
        end
        assert(~isempty(q),'Could not initialize %s',nm);
        uv=applyHomography(H,allXY(q,:));
        XYZ0(j,:)=[-65,2704-uv(1),uv(2)];
    elseif startsWith(nm,"F")
        im=5; q=find(allNames==nm & allImages==im,1);
        H=frontH;
        if isempty(q)
            im=4; q=find(allNames==nm & allImages==im,1);
            H=faceHomography(["F1" "F6" "F7" "F12"],4,allNames,allImages,allXY, ...
                [0 0;2478 0;2478 -2213;0 -2213]);
        end
        assert(~isempty(q),'Could not initialize %s',nm);
        uv=applyHomography(H,allXY(q,:));
        XYZ0(j,:)=[uv(1),0,uv(2)];
    elseif startsWith(nm,"R")
        q=find(allNames==nm & allImages==8,1);
        assert(~isempty(q),'Right-wall point %s missing in Image 8',nm);
        uv=applyHomography(rightH,allXY(q,:));
        XYZ0(j,:)=[2543,uv(1),uv(2)];
    else
        error('Unknown target %s: no face initialization rule',nm);
    end
end
assert(all(isfinite(XYZ0(:))),'Invalid XYZ initialization');

idx = startsWith(names,"R");

Rcoordinates = table( ...
    names(idx), ...
    XYZ0(idx,1), ...
    XYZ0(idx,2), ...
    XYZ0(idx,3), ...
    'VariableNames', {'Target','X_mm','Y_mm','Z_mm'});

disp(Rcoordinates);

figure;
scatter3(XYZ0(idx,1),XYZ0(idx,2),XYZ0(idx,3),60,'filled');
hold on;
text(XYZ0(idx,1),XYZ0(idx,2),XYZ0(idx,3),names(idx));
axis equal;
grid on;
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
title('Initial Right Wall Coordinates');


%% INPUT EXT: [image, camera, Xc Yc Zc omega phi kappa]
[found,rows]=ismember(activeImages,notes(:,1));assert(all(found));
chosenNotes=notes(rows,:);
ext=[chosenNotes(:,1),99*ones(numImages,1),chosenNotes(:,2:4), ...
     chosenNotes(:,7),chosenNotes(:,5),chosenNotes(:,6)];

%% GENERATE COURSE FILES
photoTable=[numericID(pointIdx),imageID,xy];
tieTable=[numericID,XYZ0];
writeRows([prefix '.pho'],photoTable,'%d %d %.4f %.4f\n');
writeRows([prefix '.tie'],tieTable,'%d %.6f %.6f %.6f\n');
writeRows([prefix '.ext'],ext,'%d %d %.6f %.6f %.6f %.6f %.6f %.6f\n');
writeRows([prefix '.int'],[99,IOP,-1],'%d %.6f %.6f %.6f %d\n');
% Include BOTH measured distances in the same nine-image bundle.
% Resolve IDs from this run's point-name map; never reuse IDs from other runs.
[tf,scaleIdx]=ismember([scaleFrom;scaleTo],names);
assert(all(tf),'F1 or F6 was removed by the minimum-image filter');
[tf2,secondIdx]=ismember([secondFrom;secondTo],names);
assert(all(tf2),'R4 or R1 was removed by the minimum-image filter');
D0=norm(XYZ0(scaleIdx(2),:)-XYZ0(scaleIdx(1),:));
D1=norm(XYZ0(secondIdx(2),:)-XYZ0(secondIdx(1),:));
distanceTable = [
    numericID(scaleIdx(1)),  numericID(scaleIdx(2)),  measuredDistance_mm;
    numericID(secondIdx(1)), numericID(secondIdx(2)), secondDistance_mm
];
writeRows([prefix '.dis'],distanceTable,'%d %d %.6f\n');
% CFG: param_obs_mode=1 (first image EOP observations); iop_mode=0 (fixed)
cfg=[1;0;0.001;0.0001;0.01;sd_xy;sd_dist;sd_firstEOP_XYZ;sd_firstEOP_wpk];
writeRows([prefix '.cfg'],cfg,'%g\n');
if isfile([prefix '.con']),delete([prefix '.con']);end

nUnknown=6*numImages+3*numPoints;
nScalar=2*numPhoto+6+size(distanceTable,1);
fprintf('\nSELECTED-IMAGE NETWORK: %d images, %d points, %d photo measurements\n', ...
    numImages,numPoints,numPhoto);
fprintf('Unknowns %d, scalar observations %d, nominal DOF %d\n',nUnknown,nScalar,nScalar-nUnknown);
fprintf('Initial F1-F6 distance: %.3f mm; tape distance: %.3f mm\n',D0,measuredDistance_mm);
fprintf('Initial R4-R1 distance: %.3f mm; tape distance: %.3f mm\n',D1,secondDistance_mm);
fprintf('Selected IDs: %s; first image %d defines EOP datum.\n', mat2str(activeImages'),activeImages(1));
fprintf('Field-note reversed phi; fixed calibrated IOP.\n');
fprintf('Interactive snooping: remove one target/image or whole target. Files: %s.*\n',prefix);
if nScalar-nUnknown<=0
    warning('Nominal degrees of freedom are nonpositive; adjustment is unlikely to be identifiable.');
end

%% RUN COURSE BUNDLE ADJUSTMENT + MANUAL DATA SNOOPING
% Requires full_matlab_bundle.p. The program prints its own convergence
% status; do not accept a returned solution if it says "did not converge".
assert(exist('full_matlab_bundle','file')==2 || exist('full_matlab_bundle','file')==6, ...
    'Add full_matlab_bundle.p to MATLAB path');

% A choice can remove one image measurement OR an entire target.
% Existing numeric IDs remain unchanged for all retained targets.
excludedObservations = strings(0,2);
% Multiplicative variance scale relative to original observation variances.
% At pass 1: cumulativeVarianceScale=1. After each valid pass:
% cumulativeVarianceScale <- cumulativeVarianceScale*(v''P_current*v/dof).
% For the next bundle, scale ALL four CFG standard deviations by sqrt(scale).
cumulativeVarianceScale = 1;
scaleHistory=zeros(0,7);  % pass, n, u, dof, local variance factor, cumulative factor, RMS
iteration = 0;
while true
    iteration = iteration+1;
    fprintf('\n====== ADJUSTMENT / DATA SNOOPING PASS %d ======\n',iteration);
    writeRows([prefix '.pho'],photoTable,'%d %d %.4f %.4f\n');
    writeRows([prefix '.tie'],tieTable,'%d %.6f %.6f %.6f\n');
    cfg(6:9)=[sd_xy;sd_dist;sd_firstEOP_XYZ;sd_firstEOP_wpk]*sqrt(cumulativeVarianceScale);
    writeRows([prefix '.cfg'],cfg,'%g\n');
    fprintf('Current cumulative variance scale = %.8f, SD multiplier = %.8f\n', ...
        cumulativeVarianceScale,sqrt(cumulativeVarianceScale));
    [eop,iop,con,tie,Cx,vr_img,vr_eop,vr_dis,vr_con] = full_matlab_bundle(prefix);

    % vr_img(:,1) is scalar residual (x,y,x,y,...); column 2 is redundancy.
    v = vr_img(:,1);
    assert(numel(v)==2*size(photoTable,1), ...
        'Image-residual count does not match PHO; cannot label residuals safely.');
    residualXY = reshape(v,2,[])';
    rmsXY=sqrt(mean(v.^2));
    fprintf('Returned image coordinate RMS = %.4f px\n',rmsXY);

    % Six weighted EOP observations, plus measured distances.
    % Number of unknowns is taken from Cx, not inferred from target count.
    numScalar = numel(v)+numel(vr_eop(:,1))+numel(vr_dis(:,1))+numel(vr_con(:,1));
    nUnknown = size(Cx,1);
    dof = numScalar-nUnknown;
    assert(dof>0,'Nonpositive degrees of freedom.');
    vPv_img=sum((v/cfg(6)).^2);
    vPv_eop=0;
    if ~isempty(vr_eop)
        % EOP order: X,Y,Z,omega,phi,kappa per camera (usually 6 for first).
        nE=numel(vr_eop(:,1));
        sigE=repmat([cfg(8);cfg(8);cfg(8);cfg(9);cfg(9);cfg(9)],ceil(nE/6),1);
        vPv_eop=sum((vr_eop(:,1)./sigE(1:nE)).^2);
    end
    vPv_dist=sum((vr_dis(:,1)/cfg(7)).^2);
    vPv_con=0;  % No .con file in this parameter-observation datum model.
    if ~isempty(vr_con)
        warning('Control residuals exist unexpectedly; include their actual SDs in vPv_con.');
    end
    vPv=vPv_img+vPv_eop+vPv_dist+vPv_con;
    sigma0sq=vPv/dof;  % local factor, using CURRENT scaled observation weights
    assert(isfinite(sigma0sq) && sigma0sq>0,'Invalid variance estimate; cannot rescale.');
    sigma0=sqrt(sigma0sq);
    % This becomes the cumulative factor of the original (unscaled) variances.
    updatedCumulativeVarianceScale=cumulativeVarianceScale*sigma0sq;
    assert(isfinite(updatedCumulativeVarianceScale) && updatedCumulativeVarianceScale>0, ...
        'Invalid cumulative scale.');
    fprintf('Cumulative variance factor (relative to original SDs) = %.8f\n',updatedCumulativeVarianceScale);
    fprintf('Cumulative SD multiplier = %.8f\n',sqrt(updatedCumulativeVarianceScale));
    scaleHistory(end+1,:)=[iteration numScalar nUnknown dof sigma0sq ...
        updatedCumulativeVarianceScale rmsXY]; %#ok<AGROW>
    writetable(array2table(scaleHistory,'VariableNames', ...
        {'Pass','Observations','Unknowns','DOF','LocalVarianceFactor', ...
         'CumulativeVarianceFactor','ImageRMS_px'}),[prefix '_variance_scale_history.csv']);
    fprintf('n=%d, u=%d, redundancy=%d, vPv=%.6f\n',numScalar,nUnknown,dof,vPv);
    fprintf('A posteriori variance factor = %.6f; sigma0 = %.6f\n',sigma0sq,sigma0);

    redundancies=vr_img(:,2);
    good=redundancies>0 & isfinite(redundancies);
    w=nan(size(v));
    % Re-standardize using the UPDATED cumulative estimated observation SD.
    % The same cumulative factor is applied to image, distance and EOP prior variances.
    w(good)=v(good)./(sd_xy*sqrt(updatedCumulativeVarianceScale).*sqrt(redundancies(good)));
    % This is a posteriori-scaled snooping; scores change with scale.
    wxy=reshape(w,2,[])';
    pointIDs=photoTable(:,1);
    [tf,pi]=ismember(pointIDs,numericID);
    assert(all(tf),'PHO IDs not found in target mapping');
    obsNames=names(pi);
    imgIDs=photoTable(:,2);
    rxy=reshape(redundancies,2,[])';
    maxAbsW=max(abs(wxy),[],2);
    maxAbsW(all(~isfinite(wxy),2))=NaN;
    residualTable=table(obsNames,pointIDs,imgIDs,photoTable(:,3),photoTable(:,4), ...
        residualXY(:,1),residualXY(:,2),rxy(:,1),rxy(:,2),wxy(:,1),wxy(:,2),maxAbsW, ...
        'VariableNames',{'Target','BundlePointID','Image','x_px','y_px', ...
        'vx_px','vy_px','rx','ry','wx','wy','MaxAbsW'});
    residualTable=sortrows(residualTable,'MaxAbsW','descend');
    crit=3;   % screening threshold; investigate before excluding
    suspects=residualTable(residualTable.MaxAbsW>crit,:);
    fprintf('Potential outlying target/image pairs (|w| > %.1f): %d\n',crit,height(suspects));
    disp(suspects(1:min(height(suspects),25),:));

    % Save the current iteration before potentially removing an observation.
    passprefix=sprintf('%s_pass%02d',prefix,iteration);
    writetable(residualTable,[passprefix '_residuals.csv']);
    save([passprefix '_results.mat'],'eop','iop','con','tie','Cx', ...
        'vr_img','vr_eop','vr_dis','vr_con','ext','IOP','cfg','names','numericID', ...
        'XYZ0','photoTable','tieTable','rmsXY','sigma0sq','sigma0','dof','vPv', ...
        'residualTable','excludedObservations','scaleHistory', ...
        'cumulativeVarianceScale','updatedCumulativeVarianceScale');

    % Prepare next pass with new weights only when a removal is confirmed.
    % Do not repeatedly multiply the scale without a new adjustment.
    % Remove a measured image point (one image) or an entire object point.
    % The original 3-image rule applies ONLY to initial selection. After
    % snooping, 2-image intersections are allowed, provided network remains
    % connected. A one-image tie point is always removed completely.
    fprintf('\nSNOOPING MENU:\n');
    fprintf('  F3 4    = remove target F3 measurement in Image 4\n');
    fprintf('  F3 all  = remove target F3 from ALL images\n');
    fprintf('  Enter   = finish and plot current adjustment\n');
    command=strtrim(string(input('Remove observation (target image, or target all): ','s')));
    if strlength(command)==0, break; end
    fields=split(command);
    if numel(fields)==1
        target=fields(1);
        fprintf('For one measurement enter e.g. "%s 4"; for entire point "%s all".\n',target,target);
        action=strtrim(string(input('Image number or ALL: ','s')));
    elseif numel(fields)==2
        target=fields(1); action=fields(2);
    else
        fprintf('Invalid entry. Use target + image number, or target + ALL.\n');
        continue;
    end
    target=upper(target);
    if ~any(names==target)
        fprintf('Target %s is not in the current network.\n',target);
        continue;
    end
    dropPoint=strcmpi(action,'all');
    if ~dropPoint
        im=str2double(action);
        if ~isscalar(im)||~isfinite(im)||~ismember(im,selectedImages)
            fprintf('Invalid image number: %s\n',action);
            continue;
        end
        ix=find(obsNames==target & imgIDs==im);
        if numel(ix)~=1
            fprintf('No measurement for %s in Image %d.\n',target,im);
            continue;
        end
        remaining=sum(obsNames==target)-1;
        if remaining<2
            fprintf('%s would have only %d image(s). Remove entire target instead.\n',target,remaining);
            dropPoint=true;
        end
    end
    if dropPoint && any(target==[scaleFrom,scaleTo,secondFrom,secondTo])
        fprintf('Cannot remove all of %s: it is an endpoint of a bundle distance.\n',target);
        continue;
    end
    if dropPoint
        removeMask=obsNames==target;
        fprintf('Removing ALL %d observations of %s and its XYZ unknown.\n',nnz(removeMask),target);
    else
        removeMask=obsNames==target & imgIDs==im;
        fprintf('Removing %s ONLY from Image %d (%d images will remain).\n',target,im,remaining);
    end
    trial=photoTable(~removeMask,:);
    if ~all(ismember(selectedImages,unique(trial(:,2))))
        fprintf('NOT removed: some image would have no observations.\n');
        continue;
    end
    % Stronger condition: keep graph of image nodes linked by common targets.
    if ~imagesConnected(trial,selectedImages)
        fprintf('NOT removed: this would disconnect the image network.\n');
        continue;
    end
    yes=strtrim(string(input('Confirm removal and rerun bundle? [y/N]: ','s')));
    if ~strcmpi(yes,'y')
        fprintf('No changes made.\n'); continue;
    end
    % Carry the a posteriori scale into the NEXT bundle's CFG weights.
    cumulativeVarianceScale=updatedCumulativeVarianceScale;
    photoTable=trial;
    if dropPoint
        pj=find(names==target,1);
        numericID(pj)=[]; names(pj)=[]; XYZ0(pj,:)=[];
        tieTable(pj,:)=[]; numPoints=numel(names);
        excludedObservations(end+1,:)=[target,"ALL"]; %#ok<AGROW>
    else
        excludedObservations(end+1,:)=[target,string(im)]; %#ok<AGROW>
    end
    writetable(table(names,numericID, ...
        'VariableNames',{'OriginalPointID','BundlePointID'}), ...
        [prefix '_target_ID_map.csv']);
    fprintf('Removed successfully; running adjustment again.\n');
end
% Last adjustment used cumulativeVarianceScale; its posterior update is
% updatedCumulativeVarianceScale, useful for FINAL reported precision.
finalCumulativeVarianceScale=updatedCumulativeVarianceScale;

% Note: do not draw statistical conclusions if the course routine reported
% "Adjustment did not converge". The P-code does not provide a convergence
% flag as a returned argument.

%% FINAL TABLES AND 3D PLOT
writetable(residualTable,[prefix '_final_residuals.csv']);
if ~isempty(excludedObservations)
    writetable(array2table(excludedObservations, ...
        'VariableNames',{'Target','Image'}),[prefix '_excluded_observations.csv']);
end

% The course uses six EOP unknowns per image, then 3 XYZ per tie point.
numEOP=size(eop.XYZc,1);
numTIE=size(tie.XYZ,1);
expectedCovSize=6*numEOP+3*numTIE;
canPlotCov=isnumeric(Cx) && isequal(size(Cx),[expectedCovSize expectedCovSize]);
if ~canPlotCov
    warning('Cx has size %dx%d; expected %d. Cannot reliably index error ellipses.', ...
        size(Cx,1),size(Cx,2),expectedCovSize);
end
assert(numTIE==numPoints,'Tie-point count changed; check output target ordering');
XYZ=tie.XYZ;
adjustedTable=table(names,XYZ(:,1),XYZ(:,2),XYZ(:,3), ...
    'VariableNames',{'Target','X_mm','Y_mm','Z_mm'});
writetable(adjustedTable,[prefix '_adjusted_points.csv']);

figure('Name','Adjusted 3D Network');hold on;grid on;axis equal;view(3);
faceColors=[0.15 .45 .85; .12 .65 .34; .91 .43 .15];
for k=1:3
    pfx=["L","F","R"];ix=startsWith(names,pfx(k));
    if any(ix)
        scatter3(XYZ(ix,1),XYZ(ix,2),XYZ(ix,3),45,faceColors(k,:), 'filled');
        text(XYZ(ix,1),XYZ(ix,2),XYZ(ix,3),cellstr(names(ix)),'FontSize',8);
    end
end
plot3(eop.XYZc(:,1),eop.XYZc(:,2),eop.XYZc(:,3),'k--^', ...
    'MarkerFaceColor','k','LineWidth',1);
for j=1:numEOP
    text(eop.XYZc(j,1),eop.XYZc(j,2),eop.XYZc(j,3), ...
        sprintf('  Camera %d',activeImages(j)),'FontSize',9);
end
xlabel('X (mm)');ylabel('Y (mm)');zlabel('Z (mm)');
title('Adjusted object points and cameras');
legend({'Left','Front','Right','Cameras'},'Location','best');
saveas(gcf,[prefix '_3D_network.png']);

%% IMAGE RESIDUAL VECTORS (MAGNIFIED FOR VISIBILITY)
rawResiduals=reshape(vr_img(:,1),2,[])';
figure('Name','Image Residual Vectors');
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
resScale=30;
for im=1:9
    nexttile;hold on;grid on;
    ok=photoTable(:,2)==im;
    scatter(photoTable(ok,3),photoTable(ok,4),8,'k','filled');
    quiver(photoTable(ok,3),photoTable(ok,4),resScale*rawResiduals(ok,1), ...
        resScale*rawResiduals(ok,2),0,'r');
    title(sprintf('Image %d',im));xlim([0 imageSize(1)]);ylim([0 imageSize(2)]);
    set(gca,'YDir','reverse');axis equal;
end
sgtitle(sprintf('Image coordinate residual vectors (%dx)',resScale));
saveas(gcf,[prefix '_residual_vectors.png']);

figure('Name','Standardized Residuals');
plot(w,'.');hold on;yline(crit,'r--');yline(-crit,'r--');grid on;
xlabel('Scalar image observation index');ylabel('Standardized residual w');
title('Image residual screening with cumulative a posteriori variance scaling');
saveas(gcf,[prefix '_standardized_residuals.png']);

%% XY ERROR ELLIPSES AND 3D CONFIDENCE ELLIPSOIDS
% Cx interpretation depends on full_matlab_bundle.p. It may be covariance
% (already scaled by the CFG variances) or a unit-weight cofactor matrix.
% Do NOT multiply it by finalCumulativeVarianceScale without verifying its
% definition. This script leaves Cx unchanged to avoid double-scaling.
% Confidence factor sqrt(chi2inv(0.95,2)) = 2.4477 for 2D;
% sqrt(chi2inv(0.95,3)) = 2.7955 for 3D.
if canPlotCov
    stdAll=sqrt(max(diag(Cx),0));
    EOPsd=reshape(stdAll(1:6*numEOP),6,[])';
    eopPrecision=array2table([activeImages(:),eop.XYZc,eop.wpk,EOPsd], ...
       'VariableNames',{'Image','X_mm','Y_mm','Z_mm','Omega_deg','Phi_deg','Kappa_deg', ...
       'sX_mm','sY_mm','sZ_mm','sOmega_deg','sPhi_deg','sKappa_deg'});
    writetable(eopPrecision,[prefix '_EOP_precision.csv']);
    XYZsd=reshape(stdAll(6*numEOP+1:end),3,[])';
    pointPrecision=table(names,XYZ(:,1),XYZ(:,2),XYZ(:,3), ...
       XYZsd(:,1),XYZsd(:,2),XYZsd(:,3), ...
       'VariableNames',{'Target','X_mm','Y_mm','Z_mm','sX_mm','sY_mm','sZ_mm'});
    writetable(pointPrecision,[prefix '_point_precision.csv']);
    ellipseScale=1;  % set >1 for display only; confidence level unchanged
    factor2=2.44774683068; factor3=2.79548348292;
    figure('Name','95% XY Error Ellipses');hold on;grid on;axis equal;
    scatter(XYZ(:,1),XYZ(:,2),15,'k','filled');
    for j=1:numTIE
        k=6*numEOP+(j-1)*3+(1:3);
        Ci=(Cx(k,k)+Cx(k,k)')/2;
        if all(isfinite(Ci(:))) && min(eig(Ci))>=-1e-9
            drawEllipse2D(Ci(1:2,1:2),XYZ(j,1:2),factor2*ellipseScale);
        end
        text(XYZ(j,1),XYZ(j,2),[' ' char(names(j))],'FontSize',8);
    end
    xlabel('X (mm)');ylabel('Y (mm)');
    title(sprintf('Point XY 95%% confidence ellipses (display scale %gx)',ellipseScale));
    saveas(gcf,[prefix '_XY_error_ellipses.png']);

    figure('Name','95% EOP Confidence Ellipsoids');hold on;grid on;axis equal;view(3);
    for j=1:numEOP
        k=(j-1)*6+(1:3);
        Ci=(Cx(k,k)+Cx(k,k)')/2;
        if all(isfinite(Ci(:))) && min(eig(Ci))>=-1e-9
            drawEllipsoid3D(Ci,eop.XYZc(j,:),factor3*5);
        end
        plot3(eop.XYZc(j,1),eop.XYZc(j,2),eop.XYZc(j,3),'k.','MarkerSize',15);
        text(eop.XYZc(j,1),eop.XYZc(j,2),eop.XYZc(j,3), ...
            sprintf(' Image %d',activeImages(j)));
    end
    xlabel('X (mm)');ylabel('Y (mm)');zlabel('Z (mm)');
    title('EOP position 95% ellipsoids (5x visualization magnification)');
    saveas(gcf,[prefix '_EOP_error_ellipsoids.png']);
end
fprintf('\nAnalysis and figures saved with prefix: %s\n',prefix);
save([prefix '_results.mat'],'eop','iop','con','tie','Cx','vr_img', ...
    'vr_eop','vr_dis','vr_con','ext','IOP','cfg','names','numericID','XYZ0', ...
    'photoTable','tieTable','rmsXY','sigma0sq','sigma0','dof','vPv', ...
    'residualTable','excludedObservations','scaleHistory', ...
    'finalCumulativeVarianceScale');

%% LOCAL FUNCTIONS
function writeRows(path,A,fmt)
fid=fopen(path,'w');assert(fid>=0,'Cannot write %s',path);
c=onCleanup(@()fclose(fid));
for j=1:size(A,1),fprintf(fid,fmt,A(j,:));end
end
function nm=pointNameFull(file)
persistent cachePath cacheNames cacheImgs cacheXY
if isempty(cachePath) || ~strcmp(cachePath,file)
 f=fopen(file,'r');assert(f>=0,'Cannot read %s',file);
 v=textscan(f,'%s %f %f %f');fclose(f);
 cacheNames=string(v{1});cacheImgs=v{2};cacheXY=[v{3},v{4}];cachePath=file;
end
nm=cacheNames;
end
function im=imageIDsFull(file)
% The persistent data in pointNameFull is purposely not shared; small PHO is cheap.
f=fopen(file,'r');v=textscan(f,'%s %f %f %f');fclose(f);im=v{2};
end
function xy=coordinatesFull(file)
f=fopen(file,'r');v=textscan(f,'%s %f %f %f');fclose(f);xy=[v{3},v{4}];
end
function H=faceHomography(corners,img,allNames,allImgs,allXY,uv)
A=zeros(8,8);b=zeros(8,1);
for t=1:4
 q=find(allNames==corners(t) & allImgs==img,1);
 assert(~isempty(q),'Corner %s missing in Image %d',corners(t),img);
 x=allXY(q,1);y=allXY(q,2);u=uv(t,1);v=uv(t,2);
 A(2*t-1,:)=[x y 1 0 0 0 -u*x -u*y];
 A(2*t,:)=[0 0 0 x y 1 -v*x -v*y];
 b(2*t-1:2*t)=[u;v];
end
assert(rcond(A)>1e-14,'Degenerate homography in Image %d',img);
h=A\b;
H=[h(1) h(2) h(3);h(4) h(5) h(6);h(7) h(8) 1];
end
function uv=applyHomography(H,xy)
p=H*[xy(:);1];assert(abs(p(3))>1e-12,'Invalid homography');uv=(p(1:2)/p(3))';
end

function drawEllipse2D(C,center,factor)
[V,D]=eig((C+C')/2);eigv=max(diag(D),0);
theta=linspace(0,2*pi,120);
xy=(V*diag(factor*sqrt(eigv))*[cos(theta);sin(theta)])';
plot(center(1)+xy(:,1),center(2)+xy(:,2),'b-','LineWidth',.8);
end

function drawEllipsoid3D(C,center,factor)
[V,D]=eig((C+C')/2);axeslength=factor*sqrt(max(diag(D),0));
[x,y,z]=sphere(15);p=[x(:) y(:) z(:)]*diag(axeslength)*V';
X=reshape(p(:,1),size(x))+center(1);
Y=reshape(p(:,2),size(y))+center(2);
Z=reshape(p(:,3),size(z))+center(3);
surf(X,Y,Z,'FaceAlpha',.2,'EdgeColor','none');
end

function connected=imagesConnected(P,imageList)
% An image graph edge exists when both images observe a retained target.
imageList=imageList(:); n=numel(imageList);
visited=false(n,1);visited(1)=true;
changed=true;
while changed
    changed=false;
    for a=1:n
        if ~visited(a),continue;end
        pts=P(P(:,2)==imageList(a),1);
        for b=1:n
            if ~visited(b) && any(ismember(pts,P(P(:,2)==imageList(b),1)))
                visited(b)=true;changed=true;
            end
        end
    end
end
connected=all(visited);
end
