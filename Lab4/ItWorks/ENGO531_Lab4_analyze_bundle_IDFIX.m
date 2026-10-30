%% ENGO 531 Lab 4 - Analyze converged 9-image bundle adjustment
% Run separately AFTER full_matlab_bundle reports convergence.
% Uses residuals, redundancy numbers, and covariance from saved results.
clear; clc; close all;
[fn,folder] = uigetfile('*_results.mat','Select converged nine-image results');
if isequal(fn,0), return; end
S = load(fullfile(folder,fn));
required = {'eop','tie','Cx','vr_img','vr_eop','vr_dis'};
missing = required(~isfield(S,required));
assert(isempty(missing), 'Results MAT is missing: %s',strjoin(missing,', '));
% The nine-image standalone script saves its .cfg separately, not in .mat.
if ~isfield(S,'cfg') || isempty(S.cfg)
    cfgPath = fullfile(folder,strrep(fn,'_results.mat','.cfg'));
    if ~isfile(cfgPath)
        [cfn,cfd] = uigetfile('*.cfg','Select matching nine-image .cfg file');
        assert(~isequal(cfn,0),'The matching .cfg file is required to obtain observation weights');
        cfgPath = fullfile(cfd,cfn);
    end
    S.cfg = readmatrix(cfgPath,'FileType','text');
end
S.cfg = S.cfg(:);
assert(numel(S.cfg)>=9 && all(isfinite(S.cfg(1:9))), ...
    'Expected nine numeric settings in matching .cfg file');
assert(size(S.vr_img,2)>=2,'Expected vr_img columns [residual, redundancy]');

nImages = numel(S.eop.img);
nTie = numel(S.tie.pt);
nImageCoordinates = size(S.vr_img,1); % x and y each count separately
nEopObs = size(S.vr_eop,1);
nDistanceObs = size(S.vr_dis,1);
n = nImageCoordinates+nEopObs+nDistanceObs;
u = 6*nImages+3*nTie;
r = n-u;
assert(r>0,'Nonpositive degrees of freedom');

% cfg: [parameter observation mode, IOP mode, convergence settings..., 
%       sigma_image_px, sigma_distance_mm, sigma_firstEOP_XYZ_mm,
%       sigma_firstEOP_wpk_deg].
sigmaImage = S.cfg(6);
sigmaDistance = S.cfg(7);
sigmaEopXYZ = S.cfg(8);
sigmaEopAngles = S.cfg(9);

% First column: residual. Second column: redundancy number.
vImage = S.vr_img(:,1);
vEop = S.vr_eop(:,1);
vDistance = S.vr_dis(:,1);
assert(nEopObs==0 || nEopObs==6, ...
    'EOP weighting below assumes six observations for first image');

imageContribution = sum((vImage/sigmaImage).^2);
distanceContribution = sum((vDistance/sigmaDistance).^2);
eopContribution = 0;
if nEopObs == 6
    % Course parameter observations correspond to XYZ and omega/phi/kappa.
    eopContribution = sum((vEop(1:3)/sigmaEopXYZ).^2) + ...
                      sum((vEop(4:6)/sigmaEopAngles).^2);
end
vTPv = imageContribution + distanceContribution + eopContribution;
varianceFactor = vTPv/r;
sigma0hat = sqrt(varianceFactor);
redundancySum = sum(S.vr_img(:,2))+sum(S.vr_eop(:,2))+sum(S.vr_dis(:,2));

fprintf('\n==== NETWORK AND A POSTERIORI VARIANCE FACTOR ====\n');
fprintf('Images: %d; tie points: %d\n',nImages,nTie);
fprintf('Image coordinates: %d; EOP observations: %d; distances: %d\n', ...
    nImageCoordinates,nEopObs,nDistanceObs);
fprintf('n = %d, u = %d, r = n-u = %d\n',n,u,r);
fprintf('Sum of redundancy numbers = %.6f (should approximately equal r)\n',redundancySum);
fprintf('v''Pv image = %.6f, EOP = %.6f, distance = %.6f\n', ...
    imageContribution,eopContribution,distanceContribution);
fprintf('v''Pv = %.6f\n',vTPv);
fprintf('A posteriori variance factor sigma0hat^2 = %.6f\n',varianceFactor);
fprintf('A posteriori unit-weight SD sigma0hat = %.6f\n',sigma0hat);
fprintf('Raw image-coordinate RMS = %.6f px\n',sqrt(mean(vImage.^2)));

%% Image residuals and standardized residuals
% Assumes image coordinates have equal a priori precision and
% the redundancy numbers correspond to their individual scalar observations.
redundancy = S.vr_img(:,2);
valid = redundancy>1e-10;
w = NaN(size(vImage));
w(valid) = vImage(valid)./(sigmaImage*sigma0hat.*sqrt(redundancy(valid)));

imageResiduals = table((1:nImageCoordinates)',vImage,redundancy,w,abs(w)>3, ...
    'VariableNames',{'ObservationCoordinateRow','Residual_px','Redundancy','StandardizedResidual','FlagAbsOver3'});
writetable(imageResiduals,fullfile(folder,'ENGO531_image_residual_analysis.csv'));
fprintf('\nStandardized image residuals |w| > 3: %d of %d coordinates\n', ...
    nnz(abs(w)>3),numel(w));

% Image observation rows are stored in the course .pho file, if present.
phoPath = fullfile(folder,replace(fn,'_results.mat','.pho'));
if ~isfile(phoPath)
    [pfn,pfd] = uigetfile('*.pho','Optional: select matching bundle .pho file (Cancel to skip)');
    if ~isequal(pfn,0),phoPath=fullfile(pfd,pfn);end
end
if isfile(phoPath)
    phot = readmatrix(phoPath,'FileType','text');
    if size(phot,1)*2 == nImageCoordinates
        imageResiduals.BundlePointID = repelem(phot(:,1),2);
        imageResiduals.ImageID = repelem(phot(:,2),2);
        imageResiduals.Axis = repmat(["x";"y"],size(phot,1),1);
        writetable(imageResiduals,fullfile(folder,'ENGO531_image_residual_analysis.csv'));
        perImage = unique(phot(:,2));
        stats = zeros(numel(perImage),5);
        for k=1:numel(perImage)
            sel = imageResiduals.ImageID==perImage(k);
            stats(k,:) = [perImage(k),sum(sel)/2,sqrt(mean(vImage(sel).^2)), ...
                         max(abs(vImage(sel))),sum(abs(w(sel))>3)];
        end
        imageStats = array2table(stats,'VariableNames', ...
            {'Image','PointObservations','RMS_perCoordinate_px','MaxAbsResidual_px','NumFlaggedCoordinates'});
        disp('PER-IMAGE RESIDUALS:');disp(imageStats);
        writetable(imageStats,fullfile(folder,'ENGO531_per_image_statistics.csv'));
    else
        warning('Matching .pho file has unexpected number of rows; skipped per-image statistics');
    end
end

%% Parameter precision from provided Cx
% As in the course Lab 1 script, parameter SD = sqrt(diagonal(Cx)).
% Verify the course's covariance convention before multiplying by sigma0hat.
Cx=S.Cx;
assert(size(Cx,1)==size(Cx,2),'Covariance matrix must be square');
assert(size(Cx,1)>=u,'Covariance is smaller than expected parameter vector');
sigmaParam=sqrt(max(0,diag(Cx)));

cameraRows = reshape(sigmaParam(1:6*nImages),6,[])';
eopPrecision=table(S.eop.img(:),cameraRows(:,1),cameraRows(:,2),cameraRows(:,3), ...
    cameraRows(:,4),cameraRows(:,5),cameraRows(:,6), ...
    'VariableNames',{'Image','StdX_mm','StdY_mm','StdZ_mm', ...
                     'StdOmega_courseUnit','StdPhi_courseUnit','StdKappa_courseUnit'});
disp('EOP PRECISION FROM Cx:');disp(eopPrecision);
writetable(eopPrecision,fullfile(folder,'ENGO531_EOP_precision.csv'));

pointRows = reshape(sigmaParam(6*nImages+(1:3*nTie)),3,[])';
tiePrecision=table(S.tie.pt(:),pointRows(:,1),pointRows(:,2),pointRows(:,3), ...
    sqrt(sum(pointRows.^2,2)), ...
    'VariableNames',{'BundlePointID','StdX_mm','StdY_mm','StdZ_mm','RSS_3D_mm'});
pointMapPath=fullfile(folder,'ENGO531_JOINT_1to9_point_map.csv');
if isfile(pointMapPath)
    pm=readtable(pointMapPath,'TextType','string');
    if all(ismember({'Target','ID'},pm.Properties.VariableNames))
        % MATLAB's bundle output can store point IDs as strings while
        % readtable imports the CSV ID column as numeric. Match as text.
        bundleKeys = strtrim(string(tiePrecision.BundlePointID));
        mapKeys = strtrim(string(pm.ID));
        [present,whichRow] = ismember(bundleKeys,mapKeys);
        labels = strings(nTie,1);
        labels(present) = string(pm.Target(whichRow(present)));
        tiePrecision.Target = labels;
        if nnz(present) ~= nTie
            warning('Matched %d of %d bundle point IDs to the point map.',nnz(present),nTie);
        end
    end
end
disp('POINT PRECISION FROM Cx:');disp(tiePrecision);
writetable(tiePrecision,fullfile(folder,'ENGO531_point_precision.csv'));

%% Plots
figure('Name','Image-coordinate residuals');
plot(vImage,'.'); grid on; hold on;
yline(3*sigmaImage,'--');yline(-3*sigmaImage,'--');
xlabel('Scalar image-coordinate observation');ylabel('Residual (pixels)');
title('Bundle adjustment image residuals');
exportgraphics(gcf,fullfile(folder,'ENGO531_image_residuals.png'),'Resolution',200);
figure('Name','Image observation redundancy');
plot(redundancy,'.');grid on;
xlabel('Scalar image-coordinate observation');ylabel('Redundancy number');
title('Image observation redundancy numbers');
exportgraphics(gcf,fullfile(folder,'ENGO531_image_redundancy.png'),'Resolution',200);

fprintf('\nSaved CSV tables and plots to:\n  %s\n',folder);
fprintf('NOTE: Coordinates and covariance require independent check distances for accuracy assessment.\n');
