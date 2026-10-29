%% ENGO 531 Lab 4 - Estimated XYZ with Original Rotation Angles
clear; clc; close all;

%% 1. Settings

imageWidth = 5712;
imageHeight = 4284;

% Fixed calibrated IOPs [xp, yp, c] (pixels)
IOP = [2834.0 2156.7 3995.9];

% Observation standard deviations
sd_xy = 3.2;         % pixels
sd_dist = 5.0;       % mm
sd_XYZ = 1.0;        % mm
sd_wpk = 0.01;       % degrees

% Scale distance
scaleFrom = "F1";
scaleTo = "F6";
measuredDistance = 2478;  % mm

% Only use connected Images 1-5
useImages = [1 2 3];

% Minimum observations per target
minImages = 3;

% F1 currently appears in only two images
allowF1Exception = true;

%% 2. Estimated camera positions and original rotation angles

% Format:
% ImageID CameraID Xc Yc Zc omega phi kappa

ext = [
    1 99 -2392  1771  -600  90 90 0;
    2 99 -2472  -107  -621  90 75 0;
    3 99 -1454 -1691  -612  90 45 0;
    4 99  -168 -2366  -522  90 25 0;
    5 99  1188 -2291  -567  90  0 0
];

%% 3. Select input file

[fn, fd] = uigetfile( ...
    {'*.PHO;*.pho','Image Observations'}, ...
    'Select calibration(2).PHO');

if isequal(fn,0)
    return;
end

phoFile = fullfile(fd,fn);

outputDir = uigetdir(fd, ...
    'Select output folder');

if isequal(outputDir,0)
    return;
end

prefix = fullfile(outputDir, ...
    'engo531_lab4_estimatedXYZ_originalAngles');

%% 4. Read PHO observations

fid = fopen(phoFile,'r');

assert(fid >= 0,'Unable to open PHO file');

C = textscan(fid,'%s %f %f %f');

fclose(fid);

pointName = string(C{1});
imageID = C{2};
xy = [C{3}, C{4}];

% Retain only Images 1-5
idx = ismember(imageID,useImages);

pointName = pointName(idx);
imageID = imageID(idx);
xy = xy(idx,:);

%% 5. Count target observations

[pointIDs,~,group] = unique(pointName,'stable');

counts = accumarray(group,1);

overlapTable = table(pointIDs,counts, ...
    'VariableNames',{'Target','ImageCount'});

disp(overlapTable);

% Retain targets appearing in at least 3 images
keepPoints = pointIDs(counts >= minImages);

% Temporary exception for F1
if allowF1Exception

    if any(pointIDs == "F1")

        if ~ismember("F1",keepPoints)

            keepPoints(end+1) = "F1";

            warning('F1 retained with only two observations');

        end
    end
end

idx = ismember(pointName,keepPoints);

pointName = pointName(idx);
imageID = imageID(idx);
xy = xy(idx,:);

[pointIDs,~,pointIndex] = unique(pointName,'stable');

nPoints = numel(pointIDs);
nImages = numel(useImages);
nObs = numel(imageID);

fprintf('\nRetained %d observations from %d targets\n', ...
    nObs,nPoints);

%% 6. Assign numeric point IDs

numericID = (1001:1000+nPoints)';

idMap = table(pointIDs,numericID, ...
    'VariableNames',{'OriginalPointID','BundlePointID'});

writetable(idMap,[prefix '_point_id_map.csv']);

%% 7. Approximate object coordinates

% Estimated initial XYZ from previous joint fit
% All coordinates in millimetres

initialNames = [
    "F1";
    "F12";
    "F3";
    "F4";
    "F6";
    "F7";
    "L10";
    "L11";
    "L2";
    "L4";
    "L5";
    "L7";
    "L9"
];

initialXYZ = [
     -50.9988    48.8390    -38.7036;
     -86.4569    51.3483  -2133.2036;
     702.9923    46.4639    -51.3504;
    1674.0878    75.8147    -70.2514;
    2425.8979   116.1343    -81.8228;
    2418.7110   149.9278  -2181.5225;
    -153.2736   591.1483  -1078.9047;
    -173.5873   735.2713  -1781.2879;
    -192.6145  1821.9242    -37.6330;
    -114.3388    16.0721    -57.5116;
    -147.9521    33.0010  -2117.0338;
    -225.3004  1839.6915  -2101.7075;
    -186.6557  1412.0734  -1071.6516
];

% Match XYZ coordinates to retained target IDs
[found,idxXYZ] = ismember(pointIDs,initialNames);

assert(all(found), ...
    'Missing initial XYZ coordinates for some targets');

XYZ0 = initialXYZ(idxXYZ,:);

%% 8. Check scale distance

[found,scaleIdx] = ismember( ...
    [scaleFrom;scaleTo],pointIDs);

assert(all(found), ...
    'Scale endpoints not retained');

initialDistance = norm( ...
    XYZ0(scaleIdx(2),:) - XYZ0(scaleIdx(1),:));

fprintf('\nInitial F1-F6 distance = %.3f mm\n', ...
    initialDistance);

fprintf('Measured F1-F6 distance = %.3f mm\n', ...
    measuredDistance);

%% 9. Write PHO file

phoOut = [
    numericID(pointIndex), ...
    imageID, ...
    xy
];

fid = fopen([prefix '.pho'],'w');

assert(fid >= 0,'Could not write PHO file');

for i = 1:size(phoOut,1)

    fprintf(fid,'%d %d %.2f %.2f\n', ...
        phoOut(i,:));

end

fclose(fid);

%% 10. Write EXT file

writematrix(ext,[prefix '.ext'], ...
    'Delimiter','tab','FileType','text');

%% 11. Write TIE file

tieOut = [numericID,XYZ0];

writematrix(tieOut,[prefix '.tie'], ...
    'Delimiter','tab','FileType','text');

%% 12. Write INT file

% Format: cameraID xp yp c y_direction
% -1 = left-handed pixel coordinate system

int = [99 IOP -1];

writematrix(int,[prefix '.int'], ...
    'Delimiter','tab','FileType','text');

%% 13. Write DIS file

dis = [
    numericID(scaleIdx(1)), ...
    numericID(scaleIdx(2)), ...
    measuredDistance
];

fid = fopen([prefix '.dis'],'w');

fprintf(fid,'%d %d %.4f\n',dis);

fclose(fid);

%% 14. Write CFG file

% Parameter observation mode = 1
% Fix IOPs = 0

cfg = [
    1;          % Six weighted EOP observations for Image 1
    0;          % IOPs fixed
    0.001;      % XYZ convergence tolerance (mm)
    0.0001;     % Angular tolerance (degrees)
    0.01;       % IOP tolerance (pixels)
    sd_xy;      % Image-coordinate SD
    sd_dist;    % Distance SD
    sd_XYZ;     % Image 1 EOP position SD
    sd_wpk      % Image 1 EOP angle SD
];

writematrix(cfg,[prefix '.cfg'], ...
    'FileType','text');

% No CON file for parameter observation mode 1
if isfile([prefix '.con'])
    delete([prefix '.con']);
end

%% 15. Network summary

nUnknowns = 6*nImages + 3*nPoints;

nObservations = 2*nObs + 6 + 1;

dof = nObservations - nUnknowns;

fprintf('\n--- Network Summary ---\n');

fprintf('Images: %d\n',nImages);
fprintf('Targets: %d\n',nPoints);
fprintf('Image points: %d\n',nObs);
fprintf('Unknowns: %d\n',nUnknowns);
fprintf('Observations: %d\n',nObservations);
fprintf('Degrees of freedom: %d\n',dof);

disp('Initial EOPs:');
disp(ext);

%% 16. Run bundle adjustment

assert(exist('full_matlab_bundle','file') == 6 || ...
       exist('full_matlab_bundle','file') == 2, ...
       'full_matlab_bundle.p must be on the MATLAB path');

[eop,iop,con,tie,Cx,vr_img,vr_eop,vr_dis,vr_con] = ...
    full_matlab_bundle(prefix);

%% 17. Results

v = vr_img(:,1);

rmsXY = sqrt(mean(v.^2));

fprintf('\n--- Bundle Adjustment Results ---\n');

fprintf('Image coordinate RMS = %.4f pixels\n',rmsXY);

disp('Adjusted EOP positions [mm]:');
disp(eop.XYZc);

disp('Adjusted EOP rotations [degrees]:');
disp(eop.wpk);

disp('Adjusted tie point coordinates [mm]:');
disp(tie.XYZ);

disp('Distance residual [mm]:');
disp(vr_dis(:,1));

%% 18. Save results

save([prefix '_results.mat'], ...
    'eop','iop','con','tie','Cx', ...
    'vr_img','vr_eop','vr_dis','vr_con', ...
    'ext','int','cfg','tieOut','phoOut', ...
    'idMap','rmsXY','dof');

%% 19. Plot residuals

figure;

plot(v,'.');

grid on;

xlabel('Residual coordinate index');
ylabel('Image residual (pixels)');

title('ENGO 531 Lab 4 - Image Coordinate Residuals');

fprintf('\nTrial completed. Check MATLAB convergence message.\n');