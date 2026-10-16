%% ENGO 531 Lab 3: prepare inputs for the COURSE full_matlab_bundle.p
% Preliminary adjustment: IOPs fixed; all observed points with known XYZ
% supplied as weighted control points. This is NOT the final 3-GCP solution.
% Make sure full_matlab_bundle.p is on the MATLAB path.
clear; clc;

%% 1. Settings: VERIFY THESE BEFORE RUNNING
imageWidthPx  = 5712;     % CHANGE to actual width of measurement JPEGs
imageHeightPx = 4284;     % CHANGE to actual height of measurement JPEGs
c0_px = 4109.10;          % preliminary principal distance in PIXELS
sd_xy = 1.0;              % assumed image-coordinate precision, pixels
sd_XYZ = 1.0;             % weighted control coordinate precision, mm
% NB: all images must have matching pixel geometry for one camera INT.

%% 2. Select your NEW Lab 3 observation file and target XYZ file
[fn,fd] = uigetfile({'*.PHO;*.pho','Image observations (*.PHO)'}, ...
    'Choose YOUR Lab 3 calibration.PHO (not the Lab 1 example)');
if isequal(fn,0), return; end
phoSource = fullfile(fd,fn);
[fn,fd] = uigetfile({'*.txt;*.con;*.csv;*.*','Target XYZ coordinates'}, ...
    'Choose END 305 target coordinates transformed.txt');
if isequal(fn,0), return; end
xyzSource = fullfile(fd,fn);

outputDir = uigetdir(pwd,'Choose folder for Lab 3 bundle input/output files');
if isequal(outputDir,0), return; end
prefix = fullfile(outputDir,'engo531_lab3_prelim');

%% 3. Read and validate observations and known object coordinates
obs = readmatrix(phoSource,'FileType','text');
xyz = readmatrix(xyzSource,'FileType','text');
assert(size(obs,2)==4,'PHO must have four columns: point_id image_id x y');
assert(size(xyz,2)==4,'XYZ must have four columns: point_id X Y Z');
assert(all(isfinite(obs(:))) && all(isfinite(xyz(:))), ...
    'Input files contain missing/non-numeric values');
assert(numel(unique(xyz(:,1))) == size(xyz,1),'Duplicate XYZ point IDs');
assert(size(unique(obs(:,1:2),'rows'),1)==size(obs,1), ...
    'Duplicate point_id/image_id pairs in PHO');
assert(all(obs(:,3)>0 & obs(:,3)<=imageWidthPx & ...
           obs(:,4)>0 & obs(:,4)<=imageHeightPx), ...
       'One or more measured pixels fall outside the specified image size');

% First three image IDs, corresponding to the EOP table provided.
% Test landscape images only
useIDs = [1 2 3 4 5 6 7 8];
obs = obs(ismember(obs(:,2),useIDs),:);
[known,idxXYZ] = ismember(obs(:,1),xyz(:,1));
if any(~known)
    fprintf('Skipping %d measurements of unknown-XYZ targets: ',sum(~known));
    fprintf('%g ',unique(obs(~known,1))); fprintf('\n');
end
obs = obs(known,:);
assert(~isempty(obs),'No observations matched known XYZ coordinates');

% The COURSE full_matlab_bundle.p requires each object/control point
% to have measurements in at least two different images (homology check).
% Keep the original source PHO untouched; filter only generated inputs.
[uniquePoints,~,group] = unique(obs(:,1));
imageCounts = accumarray(group,1);
singlePoints = uniquePoints(imageCounts < 2);
if ~isempty(singlePoints)
    fprintf('Excluding %d single-image known-XYZ points: ',numel(singlePoints));
    fprintf('%g ',singlePoints); fprintf('\n');
    obs = obs(~ismember(obs(:,1),singlePoints),:);
end
assert(~isempty(obs),'No points remain observed in at least two images');
assert(all(ismember(useIDs,unique(obs(:,2)))), ...
    'At least one of the three images has no usable overlapping points');

% Include ONLY points actually observed. In this initial weighted-control
% approach, all known points have XYZ observation equations in the bundle.
usedIDs = unique(obs(:,1));
[~,idx] = ismember(usedIDs,xyz(:,1));
con = xyz(idx,:);

%% 4. Approximate EOPs: EXT order is image camera X Y Z OMEGA PHI KAPPA
% Course EXT field order is omega phi kappa (NOT phi kappa omega).
% These three rows match the user's renamed Images 1, 2, and 3.

ext = [
    1  99  -280.7872   -4149.0255   152.2944   105    0     0;
    2  99  -137.8809   -3613.9113  1227.9180    90    0   -90;
    3  99  1311.3245   -3712.2095   415.5383   105   20     0;
    4  99  2228.7159   -3592.4673  1293.5731    90   30   -90;
    5  99  2228.7159   -2747.4673   253.5731    90   60     0;
    6  99 -1641.0098   -2594.2438  1587.8896    80  -45   -90;
    7  99  -162.8183   -2404.7969   683.2079    90    0     0;
    8  99  -162.8183   -2404.7969   683.2079    90    0   -90
];

% Filter EOPs to match selected images
ext = ext(ismember(ext(:,1), useIDs), :);
disp('EOPs included in adjustment:')
disp(ext)
%% 5. Write the course-specific text files
% PHO point_id image_id x y (two decimal places)
fid=fopen([prefix '.pho'],'w'); assert(fid>0,'Cannot write PHO');
cleaner=onCleanup(@()fclose(fid));
for i=1:size(obs,1)
    fprintf(fid,'%d %d %.2f %.2f\n',round(obs(i,1)),round(obs(i,2)),obs(i,3),obs(i,4));
end
clear cleaner;

% EXT image_id camera_id Xc Yc Zc omega phi kappa
writematrix(ext,[prefix '.ext'],'Delimiter','tab','FileType','text');
% CON point_id X Y Z (all observed targets with known XYZ)
writematrix(con,[prefix '.con'],'Delimiter','tab','FileType','text');
% INT camera_id xp yp c y_dir. Top-left pixel origin, y increases downward.
% Check principal point convention against your source images.
int = [99 imageWidthPx/2 imageHeightPx/2 c0_px -1];
writematrix(int,[prefix '.int'],'Delimiter','tab','FileType','text');

% CFG: param_obs_mode, iop_mode, tolerance_XYZ, tolerance_wpk,
% tolerance_xpypc, sd_xy, sd_dist, sd_XYZ, sd_wpk
% mode=0: weighted control coordinates; IOPs held fixed. All must >0.
cfg = [0; 0; 0.001; 0.0001; 0.01; sd_xy; 1; sd_XYZ; 1];
writematrix(cfg,[prefix '.cfg'],'FileType','text');

% No TIE or DIS file needed in mode 0. Remove stale files using this prefix.
for suffix = {'.tie','.dis'}
    stale=[prefix suffix{1}];
    if isfile(stale), delete(stale); end
end

fprintf('\nPrepared %d observations, %d known XYZ points, %d images.\n', ...
    size(obs,1),size(con,1),numel(unique(obs(:,2))));
fprintf('Input prefix: %s\n',prefix);
fprintf('Principal point: (%.2f, %.2f) px; c = %.2f px\n',int(2:4));

%% 6. Run the EXISTING COURSE adjustment executable
assert(exist('full_matlab_bundle','file')==6 || ...
       exist('full_matlab_bundle','file')==2, ...
       'Place full_matlab_bundle.p on MATLAB path / current directory.');
[eop,iop,conAdjusted,tie,Cx,vr_img,vr_eop,vr_dis,vr_con] = ...
    full_matlab_bundle(prefix);

disp('Estimated EOP structure:'); disp(eop);
disp('Returned IOP structure:'); disp(iop);

%% Examine adjusted EOPs

disp('Initial EOPs:')
disp(ext)

disp('Estimated camera coordinates:')
disp(eop.XYZc)

disp('Estimated rotation angles (omega, phi, kappa):')
disp(eop.wpk)

%% Check image residuals

figure;
plot(vr_img(:,1), '.-');
xlabel('Image-coordinate observation');
ylabel('Residual (pixels)');
title('Bundle Adjustment Image Residuals');
grid on;

if isnumeric(vr_img) && size(vr_img,2)>=1 && ~isempty(vr_img)
    rmsXY = sqrt(mean(vr_img(:,1).^2));
    fprintf('Image-coordinate residual RMS = %.3f px (per coordinate)\n',rmsXY);
    figure; plot(vr_img(:,1),'.'); grid on;
    xlabel('Image-coordinate residual index'); ylabel('Residual (pixel)');
    title('Preliminary bundle adjustment: image residuals');
end
save([prefix '_results.mat'],'eop','iop','conAdjusted','tie', ...
     'Cx','vr_img','vr_eop','vr_dis','vr_con','obs','ext','int','cfg');


% Image coordinate residuals
v = vr_img(:,1);

% Residuals grouped into x and y
v_xy = reshape(v,2,[])';

% Magnitude of residual for each target observation
v_total = hypot(v_xy(:,1),v_xy(:,2));

% Largest residuals
[sortedResiduals,idx] = sort(v_total,'descend');

disp('Largest image-point residuals (pixels):')
disp(sortedResiduals(1:min(10,end)))

%% Identify observations with largest residuals

% Read the PHO file used by the adjustment
pho = readmatrix([prefix '.pho'], 'FileType', 'text');

% Extract point IDs and image IDs
pointID = pho(:,1);
imageID = pho(:,2);

% Extract image coordinate residuals
v = vr_img(:,1);

% Group x and y residuals
v_xy = reshape(v,2,[])';

% Calculate residual magnitude
v_total = hypot(v_xy(:,1),v_xy(:,2));

% Create table
residualTable = table( ...
    pointID, imageID, ...
    v_xy(:,1), v_xy(:,2), v_total, ...
    'VariableNames', ...
    {'PointID','ImageID','Vx_pix','Vy_pix','Magnitude_pix'});

% Sort largest to smallest
residualTable = sortrows(residualTable, ...
    'Magnitude_pix','descend');

% Display 15 largest residuals
disp(residualTable(1:min(15,height(residualTable)),:));