%% ENGO 531 Lab 3: minimally constrained self-calibration
% Three weighted GCPs; remaining known-XYZ targets are tie points.
% Camera IOPs (xp, yp, c) estimated using full_matlab_bundle.p.
% This script does not overwrite preliminary adjustment files.
% MINIMUM OVERLAP: 2 images per target (instead of 3).
% CAUTION: This is an experimental sensitivity run; lecture guidance
% specifies a minimum of 3 images per target for the submitted calibration.
% Stopping at |w| <= 1.96 does not prove outliers have been eliminated.
clear; clc;
close all;

%% 1. Settings
imageWidthPx = 5712;
imageHeightPx = 4284;
c0_px = 4109.10;
sd_xy = 1.0;                 % assumed image coordinate SD (pixels)
sd_XYZ = 1.0;                % assumed control-coordinate SD (mm)
useIDs = (1:8)';

% Enter exactly THREE well-distributed target IDs, e.g. [ID1 ID2 ID3].
% Leave empty to be prompted when the script runs.
controlIDs = [1123 2031 2012];

% Set true to use previously adjusted EOPs when you have a saved
% preliminary '*_results.mat' file containing eop
useAdjustedEOPs = true;

%% 2. Choose inputs
[fn,fd] = uigetfile({'*.PHO;*.pho','Image observations (*.PHO)'}, ...
    'Select merged calibration.PHO');
if isequal(fn,0), return; end
phoSource = fullfile(fd,fn);

[fn,fd] = uigetfile({'*.txt;*.con;*.csv;*.*','Target XYZ coordinates'}, ...
    'Select END 305 target coordinate file');
if isequal(fn,0), return; end
xyzSource = fullfile(fd,fn);

outputDir = uigetdir(pwd,'Select output folder for SELF-CALIBRATION files');
if isequal(outputDir,0), return; end
prefix = fullfile(outputDir,'engo531_lab3_selfcal_snooping_min2');

% Automated data-snooping settings
criticalW = 1.96;              % two-sided standardized residual threshold
maxRemovals = 100;              % safety cap (0 = initial adjustment only)
allowDropWholeTiePoint = false; % if true, 2-image tie points may be deleted entirely
% Retained points always require at least TWO distinct images.
% Control points are protected from losing their second image.


%% 3. Validate observations and XYZ
obs = readmatrix(phoSource,'FileType','text');
xyz = readmatrix(xyzSource,'FileType','text');
assert(size(obs,2)==4 && size(xyz,2)==4, ...
    'PHO and XYZ files must each contain four numeric columns');
assert(all(isfinite(obs(:))) && all(isfinite(xyz(:))), ...
    'Input contains missing or nonnumeric entries');
assert(numel(unique(xyz(:,1)))==size(xyz,1),'Repeated XYZ target IDs');
assert(size(unique(obs(:,1:2),'rows'),1)==size(obs,1), ...
    'Repeated target/image combination in PHO');
assert(all(obs(:,3)>=1 & obs(:,3)<=imageWidthPx & ...
    obs(:,4)>=1 & obs(:,4)<=imageHeightPx), ...
    'Image coordinates exceed the specified image dimensions');

obs = obs(ismember(obs(:,2),useIDs),:);


%% Identify targets needing additional image measurements

% Read the original image observations
obs_all = readmatrix(phoSource, 'FileType', 'text');

% Only consider the 8 calibration images
obs_all = obs_all(ismember(obs_all(:,2),useIDs),:);

% Keep only targets with supplied XYZ coordinates
obs_all = obs_all(ismember(obs_all(:,1),xyz(:,1)),:);

% Count the number of images observing each target
[pointIDs,~,idx] = unique(obs_all(:,1));
nImages = accumarray(idx,1);

% Select targets observed in fewer than 4 images
targetIDs = pointIDs(nImages >= 2 & nImages < 4);

% Create a table of missing observations
missingPoint = [];
missingImage = [];
currentCount = [];

for i = 1:numel(targetIDs)

    pointID = targetIDs(i);

    % Images where target is already measured
    measuredImages = obs_all(obs_all(:,1)==pointID,2);

    % Images without measurements
    missingImages = setdiff(useIDs,measuredImages);

    for j = 1:numel(missingImages)

        missingPoint(end+1,1) = pointID;
        missingImage(end+1,1) = missingImages(j);
        currentCount(end+1,1) = numel(measuredImages);

    end
end

missingTable = table(missingPoint,missingImage,currentCount, ...
    'VariableNames',{'PointID','ImageToCheck','CurrentImageCount'});

missingTable = sortrows(missingTable, ...
    {'CurrentImageCount','PointID','ImageToCheck'});

disp(missingTable);

% Save results
writetable(missingTable,'Targets_To_Check.csv');

%% Optional pre-exclusions (leave empty for fully automatic snooping)
% Previously manually excluded observations, retained here for reference:
% [1133 8; 2010 6; 2016 6; 1104 1; 2030 6]
excludeObs = zeros(0,2);
if ~isempty(excludeObs)
    obs(ismember(obs(:,1:2),excludeObs,'rows'),:) = [];
end

[known,~] = ismember(obs(:,1),xyz(:,1));
if any(~known)
    fprintf('Excluding %d observations of targets without supplied XYZ: ',sum(~known));
    fprintf('%g ',unique(obs(~known,1))); fprintf('\n');
    % Unknown targets require reasonable initial XYZ approximations in TIE.
end
obs = obs(known,:);

%% Retain only targets observed in at least 2 images

[allPts,~,grp] = unique(obs(:,1));

% Count observations for each target
counts = accumarray(grp,1);

% Identify targets observed in fewer than 2 images
removePts = allPts(counts < 2);

if ~isempty(removePts)
    fprintf('Excluding %d targets observed in fewer than 2 images: ', ...
        numel(removePts));
    fprintf('%g ',removePts);
    fprintf('\n');

    % Remove all observations belonging to those targets
    obs = obs(~ismember(obs(:,1),removePts),:);
end
assert(~isempty(obs),'No remaining multi-image observations');
assert(isequal(sort(unique(obs(:,2))),sort(useIDs)), ...
    'At least one image has no usable measurements');

usedIDs = unique(obs(:,1));
[~,ii] = ismember(usedIDs,xyz(:,1));
usedXYZ = xyz(ii,:);

%% 4. Select exactly three well-distributed control points
fprintf('\nRetained targets eligible for control selection:\n');
disp(array2table(usedXYZ,'VariableNames',{'PointID','X_mm','Y_mm','Z_mm'}));
if isempty(controlIDs)
    answer = inputdlg({'Enter exactly THREE point IDs (e.g. 1115 1133 2022):'}, ...
        'Select three weighted GCPs',[1 70]);
    if isempty(answer), return; end
    controlIDs = sscanf(strrep(answer{1},',',' '),'%f')';
end
controlIDs=controlIDs(:);
assert(numel(controlIDs)==3 && numel(unique(controlIDs))==3 && ...
       all(ismember(controlIDs,usedIDs)), ...
       'Control IDs must be three DISTINCT points retained in the PHO');

isControl = ismember(usedXYZ(:,1),controlIDs);
con = usedXYZ(isControl,:);          % point_id X Y Z : weighted GCPs
tieInit = usedXYZ(~isControl,:);     % point_id X Y Z : INITIAL guesses
assert(size(con,1)==3,'Exactly three controls required');
assert(size(tieInit,1)>0,'No tie points remain');

% A basic geometric warning, not a substitute for inspecting ZOD / network.
xyz3=con(:,2:4);
triangleArea3D=0.5*norm(cross(xyz3(2,:)-xyz3(1,:), ...
                                xyz3(3,:)-xyz3(1,:)));
assert(triangleArea3D>1,'Chosen controls are nearly collinear in XYZ');
fprintf('Selected control IDs: %s\n',num2str(con(:,1)'));
fprintf('3D control-triangle area: %.1f mm^2\n',triangleArea3D);
fprintf('%d tie points will have freely adjusted XYZ.\n',size(tieInit,1));

%% 5. Initial exterior orientation parameters
ext = [
    1 99  -280.7872 -4149.0255  152.2944 105   0    0;
    2 99  -137.8809 -3613.9113 1227.9180  90   0  -90;
    3 99  1311.3245 -3712.2095  415.5383 105  20    0;
    4 99  2228.7159 -3592.4673 1293.5731  90  30  -90;
    5 99  2228.7159 -2747.4673  253.5731  90  60    0;
    6 99 -1641.0098 -2594.2438 1587.8896  80 -45  -90;
    7 99  -162.8183 -2404.7969  683.2079  90   0    0;
    8 99  -162.8183 -2404.7969  683.2079  90   0  -90
];
ext = ext(ismember(ext(:,1),useIDs),:);

if useAdjustedEOPs
    [fn,fd] = uigetfile({'*_results.mat;*.mat','Saved previous bundle results'}, ...
        'Select your successful FIXED-IOP results MAT file (Cancel to use field EOPs)');
    if ~isequal(fn,0)
        old=load(fullfile(fd,fn),'eop');
        assert(isfield(old,'eop'),'MAT file must contain an eop structure');
        ids=str2double(string(old.eop.img));
        assert(all(isfinite(ids)),'Previous EOP image IDs could not be read');
        for k=1:size(ext,1)
            j=find(ids==ext(k,1),1);
            assert(~isempty(j),'Saved EOP missing image %d',ext(k,1));
            ext(k,3:5)=old.eop.XYZc(j,:);
            ext(k,6:8)=old.eop.wpk(j,:);
        end
        fprintf('Using adjusted EOPs from %s\n',fn);
    end
end

%% 6-7. Iterative bundle adjustment and automated data snooping
% Each image point represents TWO coordinate residuals (x and y).
% At each pass: estimate posterior unit-weight SD; calculate w-values using
% posterior-scaled standard deviations; remove at most ONE image point.
% Scale BOTH image and control SDs to preserve their relative weights.
% A 1.96 cutoff is aggressive with hundreds of simultaneous tests.

assert(exist('full_matlab_bundle','file')==6 || ...
       exist('full_matlab_bundle','file')==2, ...
       'Add full_matlab_bundle.p to the MATLAB path');

int = [99 imageWidthPx/2 imageHeightPx/2 c0_px -1];
sd_xy_run = sd_xy;
sd_XYZ_run = sd_XYZ;
snoopingLog = table('Size',[0 10], ...
    'VariableTypes',{'double','double','double','double','double', ...
                     'double','double','double','double','string'}, ...
    'VariableNames',{'Iteration','PointID','ImageID','wx','wy', ...
                     'MaxAbsW','Sigma0','SDxyBefore','ObservationsRemoved','Action'});
snoopingHistory = table('Size',[0 7], ...
    'VariableTypes',repmat({'double'},1,7), ...
    'VariableNames',{'Iteration','ImagePoints','TiePoints','DOF', ...
                     'Sigma0Sq','SDxyBefore','MaxAbsW'});
stopReason = "";

for snoopIter = 0:maxRemovals
    % Rebuild tie-point list after any dropped targets.
    usedIDs = unique(obs(:,1));
    assert(all(ismember(controlIDs,usedIDs)), ...
        'Data snooping removed one of the three datum control points');
    [~,ii] = ismember(usedIDs,xyz(:,1));
    usedXYZ = xyz(ii,:);
    con = usedXYZ(ismember(usedXYZ(:,1),controlIDs),:);
    tieInit = usedXYZ(~ismember(usedXYZ(:,1),controlIDs),:);
    assert(size(con,1)==3 && ~isempty(tieInit), ...
        'Insufficient points for weighted-three-control adjustment');

    % Use a unique prefix for EACH adjustment: history is never overwritten.
    runPrefix = sprintf('%s_iter_%02d',prefix,snoopIter);
    fid = fopen([runPrefix '.pho'],'w');
    assert(fid>=0,'Unable to write observation file');
    for k=1:size(obs,1)
        fprintf(fid,'%d %d %.2f %.2f\n',obs(k,1),obs(k,2),obs(k,3),obs(k,4));
    end
    fclose(fid);
    writematrix(ext,[runPrefix '.ext'],'Delimiter','tab','FileType','text');
    writematrix(con,[runPrefix '.con'],'Delimiter','tab','FileType','text');
    writematrix(tieInit,[runPrefix '.tie'],'Delimiter','tab','FileType','text');
    writematrix(int,[runPrefix '.int'],'Delimiter','tab','FileType','text');
    cfg = [0; 1; 0.001; 0.0001; 0.01; sd_xy_run; 1; sd_XYZ_run; 1];
    writematrix(cfg,[runPrefix '.cfg'],'FileType','text');
    if isfile([runPrefix '.dis']), delete([runPrefix '.dis']); end

    try
        [eop,iop,conAdjusted,tie,Cx,vr_img,vr_eop,vr_dis,vr_con] = ...
            full_matlab_bundle(runPrefix);
    catch ME
        error('Adjustment failed at iteration %d. Previous iteration files remain available.\n%s', ...
            snoopIter,ME.message);
    end

    assert(numel(vr_img(:,1))==2*size(obs,1), ...
        'Mismatch between image residuals and PHO observation order');
    dof_iter = size(vr_img,1)+size(vr_con,1)-size(Cx,1);
    assert(dof_iter>0,'Nonpositive degrees of freedom');
    sigma0_sq_iter = (sum(vr_img(:,1).^2)/sd_xy_run^2 + ...
                      sum(vr_con(:,1).^2)/sd_XYZ_run^2)/dof_iter;
    assert(isfinite(sigma0_sq_iter) && sigma0_sq_iter>0, ...
        'Invalid a posteriori variance factor');
    sigma0_iter = sqrt(sigma0_sq_iter);

    r_img = vr_img(:,2);
    assert(all(isfinite(r_img)) && all(r_img>0), ...
        'Invalid/nonpositive image redundancy number');
    % Posterior-standardized w: v / (sd_prior * sigma0 * sqrt(r))
    w_img = vr_img(:,1)./(sd_xy_run*sigma0_iter*sqrt(r_img));
    w_xy = reshape(w_img,2,[])';
    w_max = max(abs(w_xy),[],2);

    snoopingHistory(end+1,:) = {snoopIter,size(obs,1),size(tieInit,1), ...
        dof_iter,sigma0_sq_iter,sd_xy_run,max(w_max)};
    fprintf('Snoop %2d: %3d image points; DOF=%d; sigma0^2=%.4f; ', ...
        snoopIter,size(obs,1),dof_iter,sigma0_sq_iter);
    fprintf('posterior SDxy=%.4f px; max |w|=%.4f\n', ...
        sd_xy_run*sigma0_iter,max(w_max));

    % Keep full data for the solution from THIS pass.
    save([runPrefix '_results.mat'],'eop','iop','conAdjusted','tie', ...
        'Cx','vr_img','vr_eop','vr_dis','vr_con','obs','ext','int', ...
        'cfg','con','tieInit','w_xy','sigma0_sq_iter','snoopingLog');

    if max(w_max)<=criticalW
        stopReason = "All posterior-standardized residuals within cutoff";
        break;
    end
    if snoopIter==maxRemovals
        stopReason = "Reached maximum removal limit";
        break;
    end

    % Descending candidates: protect control points and 2-image targets.
    [~,candidateOrder] = sort(w_max,'descend');
    chosen = 0; dropWholePoint = false;
    for m=1:numel(candidateOrder)
        k=candidateOrder(m);
        if w_max(k)<=criticalW, break; end
        pt = obs(k,1);
        ptCount = sum(obs(:,1)==pt);
        if ptCount>2
            chosen=k; break;
        elseif ptCount==2 && ~ismember(pt,controlIDs) && allowDropWholeTiePoint
            chosen=k; dropWholePoint=true; break;
        end
    end
    if chosen==0
        stopReason = "Only protected 2-image targets exceed cutoff";
        break;
    end

    badPoint = obs(chosen,1);
    badImage = obs(chosen,2);
    if dropWholePoint
        removeMask = obs(:,1)==badPoint;
        action = "Dropped entire 2-image tie point";
    else
        removeMask = false(size(obs,1),1);
        removeMask(chosen)=true;
        action = "Dropped one image point (x and y)";
    end
    % Do not allow an image to become empty.
    prospectiveObs = obs(~removeMask,:);
    if ~isequal(sort(unique(prospectiveObs(:,2))),sort(useIDs))
        stopReason = "Removal would leave an image with no measurements";
        break;
    end
    fprintf('  Removing Point %d Image %d (wx=%.3f, wy=%.3f): %s\n', ...
        badPoint,badImage,w_xy(chosen,1),w_xy(chosen,2),action);
    snoopingLog(end+1,:) = {snoopIter,badPoint,badImage, ...
        w_xy(chosen,1),w_xy(chosen,2),w_max(chosen), ...
        sigma0_iter,sd_xy_run,sum(removeMask),action};
    obs = prospectiveObs;
    % Confirm every retained target remains observed in at least 2 images.
    [remainingPts,~,remainingGrp] = unique(obs(:,1));
    remainingCounts = accumarray(remainingGrp,1);
    assert(all(remainingCounts >= 2), ...
        'Snooping left a retained target in fewer than 2 images');

    % Variance-component rescaling: common multiplier of all SDs.
    % This preserves relative weights and avoids inconsistent sd_XYZ/sd_xy.
    sd_xy_run = sd_xy_run*sigma0_iter;
    sd_XYZ_run = sd_XYZ_run*sigma0_iter;
    % Use latest EOP and IOP estimates as next approximations.
    ids = str2double(string(eop.img));
    for k=1:size(ext,1)
        j=find(ids==ext(k,1),1);
        ext(k,3:5)=eop.XYZc(j,:);
        ext(k,6:8)=eop.wpk(j,:);
    end
    int(2:4)=reshape(iop.xp_yp_c,1,[]);
end

fprintf('\nAutomated snooping finished: %s\n',stopReason);
fprintf('Total removed image points: %d\n',sum(snoopingLog.ObservationsRemoved));
fprintf('Final image points: %d\n',size(obs,1));
fprintf('Final fitted posterior-scaled image SD: %.4f px\n', ...
    sd_xy_run*sigma0_iter);
% Downstream code uses the SDs and corresponding final adjustment.
sd_xy = sd_xy_run;
sd_XYZ = sd_XYZ_run;
% Final parameter/observation files are available under runPrefix.
writetable(snoopingHistory,[prefix '_snooping_history.csv']);
writetable(snoopingLog,[prefix '_snooping_removals.csv']);
save([prefix '_snooping_summary.mat'],'snoopingHistory','snoopingLog', ...
    'stopReason','criticalW','sd_xy','sd_XYZ');

fprintf('\nInitial IOPs (xp,yp,c) [px]:\n'); disp(int(2:4));
fprintf('Estimated IOPs (xp,yp,c) [px]:\n'); disp(iop.xp_yp_c);
fprintf('Adjusted camera XYZ [mm]:\n'); disp(eop.XYZc);
fprintf('Adjusted camera wpk [degrees]:\n'); disp(eop.wpk);

%% 8. Residual statistics and largest discrepancies
v=vr_img(:,1);
assert(numel(v)==2*size(obs,1), ...
    'Unexpected residual vector length: verify correspondence to PHO');
rmsXY=sqrt(mean(v.^2));
fprintf('Image coordinate residual RMS (per coordinate) = %.3f px\n',rmsXY);

% Assumes course vr_img is ordered x,y for each input PHO measurement.
V=reshape(v,2,[])';
residualTable=table(obs(:,1),obs(:,2),V(:,1),V(:,2),hypot(V(:,1),V(:,2)), ...
    'VariableNames',{'PointID','ImageID','Vx_pix','Vy_pix','Magnitude_pix'});
residualTable=sortrows(residualTable,'Magnitude_pix','descend');
disp('Largest residuals:');
disp(residualTable(1:min(15,height(residualTable)),:));

imageIDs=unique(obs(:,2));
imageRMS=zeros(numel(imageIDs),1);
for k=1:numel(imageIDs)
    ix=obs(:,2)==imageIDs(k);
    vals=V(ix,:);
    imageRMS(k)=sqrt(mean(vals(:).^2));
end
perImageRMS=table(imageIDs,imageRMS,'VariableNames',{'ImageID','RMS_pix'});
disp('Per-image residual RMS:'); disp(perImageRMS);

save([prefix '_results.mat'],'eop','iop','conAdjusted','tie','Cx', ...
    'vr_img','vr_eop','vr_dis','vr_con','obs','ext','int','cfg', ...
    'con','tieInit','residualTable','perImageRMS','rmsXY');

figure;
plot(v,'.'); grid on;
xlabel('Residual coordinate index'); ylabel('Residual (pixels)');
title('Three-control-point self-calibration image residuals');

% Assumed parameter ordering:
% 48 EOPs, 3 IOPs, 165 object coordinates

C_IOP = Cx(49:51,49:51);

% Standard deviations
sigma_IOP = sqrt(diag(C_IOP));

% Correlation matrix
R_IOP = C_IOP ./ (sigma_IOP * sigma_IOP');

disp('IOP standard deviations:')
disp(sigma_IOP)

disp('IOP correlation matrix:')
disp(R_IOP)


%% 9. A posteriori variance factor

% Number of observations
n_img = numel(vr_img(:,1));
n_con = numel(vr_con(:,1));

% Number of estimated unknowns
u = size(Cx,1);

% Degrees of freedom
dof = n_img + n_con - u;

% Weighted sum of squared residuals
vPv_img = sum(vr_img(:,1).^2) / sd_xy^2;
vPv_con = sum(vr_con(:,1).^2) / sd_XYZ^2;

% A posteriori variance factor
sigma0_sq = (vPv_img + vPv_con) / dof;

% A posteriori standard deviation of unit weight
sigma0 = sqrt(sigma0_sq);

fprintf('Degrees of freedom: %d\n',dof);
fprintf('A posteriori variance factor: %.4f\n',sigma0_sq);
fprintf('A posteriori SD of unit weight: %.4f\n',sigma0);

%% 10. Standardized residuals and outlier detection

% A priori image-coordinate standard deviation (pixels)
sigma_xy = sd_xy * sigma0; % posterior-scaled SD for final snooping

% Residuals and redundancy numbers
v_img = vr_img(:,1);
r_img = vr_img(:,2);

% Calculate standardized residuals
w_img = v_img ./ (sigma_xy .* sqrt(r_img));

% Separate x and y standardized residuals
w_xy = reshape(w_img,2,[])';

% Critical value (95% two-sided test)
crit = 1.96;

% Create outlier table
outlierTable = table( ...
    obs(:,1), obs(:,2), ...
    w_xy(:,1), w_xy(:,2), ...
    'VariableNames', ...
    {'PointID','ImageID','wx','wy'});

% Identify outliers
outlierTable.OutlierX = abs(outlierTable.wx) > crit;
outlierTable.OutlierY = abs(outlierTable.wy) > crit;

% Keep observations flagged in either coordinate
outliers = outlierTable( ...
    outlierTable.OutlierX | outlierTable.OutlierY, :);

% Sort by largest standardized residual
outliers.MaxAbsW = max(abs(outliers.wx),abs(outliers.wy));
outliers = sortrows(outliers,'MaxAbsW','descend');

disp('Potential image-coordinate outliers:')
disp(outliers)

fprintf('Number of flagged image points: %d\n', ...
    height(outliers));

%% Analyze image overlap for each target

[pointIDs,~,idx] = unique(obs(:,1));

nImages = accumarray(idx,1);

overlapTable = table(pointIDs,nImages, ...
    'VariableNames',{'PointID','NumberOfImages'});

disp('Number of targets by image overlap:')

overlapSummary = groupsummary(overlapTable,'NumberOfImages');

disp(overlapSummary)

figure;
histogram(nImages, ...
    'BinMethod','integers');

xlabel('Number of images observing each target');
ylabel('Number of targets');
title('Image overlap distribution');
grid on;


%% 11. Plot EOP confidence ellipsoids

% Settings
S = 5;             % Ellipsoid scale factor
n = 20;            % Number of ellipsoid surface divisions
confidence = sqrt(7.8147); % 95% confidence, 3D (chi-square, 3 DOF)

figure;
hold on;
grid on;
axis equal;
view(3);

for i = 1:eop.num_eop

    % Adjusted camera coordinates
    XYZ = eop.XYZc(i,:);

    % Covariance matrix for camera position (X,Y,Z)
    idx = (i-1)*6 + (1:3);
    Ci = Cx(idx,idx);

    % Plot the 95% confidence ellipsoid
    calc_plot_ellipsoid(Ci,XYZ,confidence*S,n);

    % Plot camera centre
    plot3(XYZ(1),XYZ(2),XYZ(3),'k.', ...
        'MarkerSize',18);

    % Label image
    text(XYZ(1),XYZ(2),XYZ(3), ...
        sprintf('  Image %s',string(eop.img(i))), ...
        'FontSize',10);

end

xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
title('Adjusted EOPs with 95% Confidence Ellipsoids');

colorbar;
hold off;



%% 14. Image-point residual vectors for all images

% Convert residual vector into x and y components
vxy = reshape(vr_img(:,1),2,[])';

% Visual exaggeration of residual vectors
scale = 30;

% Create figure with 8 subplots
figure;
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');

for imageID = useIDs'

    % Select observations for current image
    idx = obs(:,2) == imageID;

    % Measured image coordinates
    xy = obs(idx,3:4);

    % Corresponding residuals
    v = vxy(idx,:);

    nexttile;
    hold on;
    grid on;
    axis equal;

    % Plot image points
    plot(xy(:,1),xy(:,2),'k.','MarkerSize',8);

    % Plot residual vectors
    quiver(xy(:,1),xy(:,2), ...
           scale*v(:,1),scale*v(:,2), ...
           0,'r','LineWidth',1.2);

    % Image coordinate system
    set(gca,'YDir','reverse');

    xlim([0 imageWidthPx]);
    ylim([0 imageHeightPx]);

    xlabel('x (pixels)');
    ylabel('y (pixels)');
    title(sprintf('Image %d',imageID));

    hold off;
end

sgtitle(sprintf('Image Coordinate Residual Vectors (%dx Exaggerated)',scale));

%% Correlation between IOPs and EOPs

% Full parameter correlation matrix
sigma = sqrt(diag(Cx));
R = Cx ./ (sigma*sigma');

% First 48 parameters are EOPs
% Parameters 49:51 are IOPs
R_IOP_EOP = R(49:51,1:48);

% Maximum absolute correlation for each IOP
[maxCorr,idx] = max(abs(R_IOP_EOP),[],2);

eopNames = {'Xc','Yc','Zc','omega','phi','kappa'};
iopNames = {'xp','yp','c'};

for i = 1:3

    imageID = ceil(idx(i)/6);
    paramID = mod(idx(i)-1,6)+1;

    fprintf('%s: Image %d, %s, correlation = %.4f\n', ...
        iopNames{i},imageID,eopNames{paramID}, ...
        R_IOP_EOP(i,idx(i)));
end