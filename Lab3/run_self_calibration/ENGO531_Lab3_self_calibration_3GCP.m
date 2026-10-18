%% ENGO 531 Lab 3: nearly minimally constrained self-calibration
% Three weighted GCPs; remaining known-XYZ targets are tie points.
% Camera IOPs (xp, yp, c) estimated using full_matlab_bundle.p.
% This script does not overwrite your preliminary adjustment files.
clear; clc;

%% 1. Settings
imageWidthPx = 5712;
imageHeightPx = 4284;
c0_px = 4109.10;
sd_xy = 3.28;                 % assumed image coordinate SD (pixels)
sd_XYZ = 1.0;                % assumed control-coordinate SD (mm)
useIDs = (1:8)';

% Enter exactly THREE well-distributed target IDs, e.g. [ID1 ID2 ID3].
% Leave empty to be prompted when the script runs.
controlIDs = [1123 2031 2012];

% Set true to use previously adjusted EOPs when you have a saved
% preliminary '*_results.mat' file containing eop (recommended).
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
prefix = fullfile(outputDir,'engo531_lab3_selfcal_3gcp');

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

%% Temporary outlier exclusion

excludeObs = [
    1133 8;
    2010 6;
    2016 6;
];

idxRemove = ismember(obs(:,1:2),excludeObs,'rows');

obs(idxRemove,:) = [];


[known,~] = ismember(obs(:,1),xyz(:,1));
if any(~known)
    fprintf('Excluding %d observations of targets without supplied XYZ: ',sum(~known));
    fprintf('%g ',unique(obs(~known,1))); fprintf('\n');
    % Unknown targets require reasonable initial XYZ approximations in TIE.
end
obs = obs(known,:);

% The course program's homology check requires two or more observations
% of each retained target from distinct images.
[allPts,~,grp] = unique(obs(:,1));
counts = accumarray(grp,1);
singlePts = allPts(counts<2);
if ~isempty(singlePts)
    fprintf('Excluding %d single-image targets: ',numel(singlePts));
    fprintf('%g ',singlePts); fprintf('\n');
    obs = obs(~ismember(obs(:,1),singlePts),:);
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

%% 6. Write bundle input files
fid=fopen([prefix '.pho'],'w');
assert(fid>=0,'Unable to create PHO');
closer=onCleanup(@()fclose(fid));
for k=1:size(obs,1)
    fprintf(fid,'%d %d %.2f %.2f\n',obs(k,1),obs(k,2),obs(k,3),obs(k,4));
end
clear closer;

writematrix(ext,[prefix '.ext'],'Delimiter','tab','FileType','text');
writematrix(con,[prefix '.con'],'Delimiter','tab','FileType','text');
writematrix(tieInit,[prefix '.tie'],'Delimiter','tab','FileType','text');
int = [99 imageWidthPx/2 imageHeightPx/2 c0_px -1];
writematrix(int,[prefix '.int'],'Delimiter','tab','FileType','text');

% CFG rows: parameter observation mode; IOP mode; XYZ, angle, IOP
% convergence tolerances; image, distance, control XYZ, EOP SDs.
% Mode 0 = weighted control coordinates; mode 1 = adjust xp, yp and c.
cfg = [0; 1; 0.001; 0.0001; 0.01; sd_xy; 1; sd_XYZ; 1];
writematrix(cfg,[prefix '.cfg'],'FileType','text');
if isfile([prefix '.dis']), delete([prefix '.dis']); end

fprintf('\nPrepared %d image measurements from %d images.\n', ...
    size(obs,1),numel(unique(obs(:,2))));
fprintf('%d weighted GCPs; %d tie points; IOP estimation enabled.\n', ...
    size(con,1),size(tieInit,1));
fprintf('Input file prefix: %s\n',prefix);

%% 7. Run the course bundle adjustment
assert(exist('full_matlab_bundle','file')==6 || ...
       exist('full_matlab_bundle','file')==2, ...
       'Add full_matlab_bundle.p to the MATLAB path');
[eop,iop,conAdjusted,tie,Cx,vr_img,vr_eop,vr_dis,vr_con] = ...
    full_matlab_bundle(prefix);

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
sigma_xy = sd_xy;

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


%% 14. Image-point residual vectors

imageID = 1;

idx = obs(:,2) == imageID;

xy = obs(idx,3:4);
vxy = reshape(vr_img(:,1),2,[])';
v = vxy(idx,:);

% Visual exaggeration of residual vectors
scale = 20;

figure;
hold on;
grid on;
axis equal;

quiver(xy(:,1),xy(:,2), ...
       scale*v(:,1),scale*v(:,2), ...
       0,'r','LineWidth',1.3);

plot(xy(:,1),xy(:,2),'k.','MarkerSize',10);

set(gca,'YDir','reverse');

xlim([0 imageWidthPx]);
ylim([0 imageHeightPx]);

xlabel('x (pixels)');
ylabel('y (pixels)');
title(sprintf('Image %d: Residual Vectors (20x Exaggerated)',imageID));

hold off;