%% ENGO 531 Lab 4 - R10 single-measurement localization diagnostic
% Select UPDATED PHO. Runs the SAME input generation as selectable-image
% convergence script, using controlled one-observation replacement trials.
% Requires the course full_matlab_bundle.p on MATLAB path.
% Original PHO is NEVER changed. All units mm and pixels as course files.
clear; clc; close all;

%% EDIT THESE OPTIONS
selectedImages = 1:7;
minImagesPerPoint = 2;   % R targets in the 1:7 subset have only 2 images
% R10 appears in exactly two selected images (6 and 7). Removing only its
% Image 6 or Image 7 observation leaves a one-image tie point, which the
% network filter removes. Those two deletion tests would therefore be the
% SAME as deleting R10 everywhere and cannot diagnose which observation
% is inconsistent. Instead compare controlled replacements using other
% shared coplanar R points to predict R10's image coordinate.
% Replacements are SYNTHETIC and are NOT valid for the final adjustment.
testNames = ["baseline","without_R10", ...
             "replace_R10_image7","replace_R10_image6"];
testExclusions = {strings(0,1), "R10", strings(0,1), strings(0,1)};
testReplacementImage = [0,0,7,6];

[fn,fd] = uigetfile({'*.PHO;*.pho','PHO point measurements'}, ...
    'Select updated calibration(4).PHO');
if isequal(fn,0), return; end
phoFile=fullfile(fd,fn);
outDir=uigetdir(fd,'Select folder for diagnostic output files');
if isequal(outDir,0), return; end
assert(numel(testNames)==numel(testExclusions),'Test names and exclusions count differ');
assert(ismember(6,selectedImages) && ismember(7,selectedImages), ...
    'Include Images 6 and 7 for the right-wall homography diagnostic');

%% PLANAR CHECK: independent of camera EOP/XYZ starting estimates
% For each common R target, fit a projective mapping from Image 6 to Image 7
% using OTHER shared R targets and calculate held-out discrepancy.
% This is a diagnostic ONLY: it assumes all targets on one plane.
fid=fopen(phoFile,'r'); assert(fid>=0,'Cannot open PHO');
V=textscan(fid,'%s %f %f %f'); fclose(fid);
pn=string(V{1}); im=V{2}; pts=[V{3},V{4}];
right=unique(pn(startsWith(pn,"R") & im==6 & ismember(pn,pn(im==7))));
right=right(arrayfun(@(x) any(pn==x & im==7),right));
A=zeros(numel(right),2);B=A;
for j=1:numel(right)
    A(j,:)=pts(find(pn==right(j)&im==6,1),:);
    B(j,:)=pts(find(pn==right(j)&im==7,1),:);
end
err=NaN(numel(right),1);
for j=1:numel(right)
    keep=true(numel(right),1);keep(j)=false;
    if sum(keep)<4,continue;end
    H=fitPlanarH(A(keep,:),B(keep,:));
    v=H*[A(j,:) 1]';
    if isfinite(v(3)) && abs(v(3))>1e-12
        pred=(v(1:2)/v(3))'; err(j)=norm(pred-B(j,:));
    end
end
Htable=table(right,err,'VariableNames',{'Target','HeldOutHomographyError_px'});
Htable=sortrows(Htable,'HeldOutHomographyError_px','descend');
fprintf('\nLEFT-OUT PLANAR TRANSFER CHECK (Image 6 -> Image 7):\n');
disp(Htable);
writetable(Htable,fullfile(outDir,'image6_to_7_homography_diagnostic.csv'));

%% PREDICT R10 IN EACH IMAGE FROM ALL OTHER COMMON R TARGETS
j10=find(right=="R10",1);
assert(~isempty(j10),'R10 is not present in both Image 6 and Image 7');
other=true(numel(right),1); other(j10)=false;
H67=fitPlanarH(A(other,:),B(other,:));
H76=fitPlanarH(B(other,:),A(other,:));
q7=H67*[A(j10,:) 1]'; predicted7=(q7(1:2)/q7(3))';
q6=H76*[B(j10,:) 1]'; predicted6=(q6(1:2)/q6(3))';
assert(all(isfinite([predicted6 predicted7])),'Invalid projected R10');
fprintf('\nR10 measured Image 6: [%.3f %.3f], synthetic prediction: [%.3f %.3f], discrepancy %.3f px\n', ...
 A(j10,:),predicted6,norm(A(j10,:)-predicted6));
fprintf('R10 measured Image 7: [%.3f %.3f], synthetic prediction: [%.3f %.3f], discrepancy %.3f px\n', ...
 B(j10,:),predicted7,norm(B(j10,:)-predicted7));

%% RUN COURSE ADJUSTMENT FOR EACH DIAGNOSTIC
summary=struct('Trial',{},'Excluded',{},'RMS_px',{},'MaxCameraShift_mm',{}, ...
    'Finite',{},'Status',{},'Message',{});
for t=1:numel(testNames)
    fprintf('\n====================================================\n');
    fprintf('TEST %d/%d: %s\n',t,numel(testNames),testNames(t));
    fprintf('====================================================\n');
    try
        repImage=testReplacementImage(t);
        repXY=[];
        if repImage==6, repXY=predicted6; end
        if repImage==7, repXY=predicted7; end
        r=runOneTrial(phoFile,outDir,selectedImages,minImagesPerPoint, ...
            string(testExclusions{t}),char(testNames(t)),repImage,repXY);
    catch ME
        r=struct('Trial',testNames(t), ...
          'Excluded',strjoin(string(testExclusions{t}),','),'RMS_px',NaN, ...
          'MaxCameraShift_mm',NaN,'Finite',false,'Status',"ERROR", ...
          'Message',string(ME.message));
        warning('Trial %s failed: %s',testNames(t),ME.message);
    end
    summary(end+1)=r; %#ok<SAGROW>
end
T=struct2table(summary);
fprintf('\n===== SUMMARY (check MATLAB convergence messages above) =====\n');
disp(T);
writetable(T,fullfile(outDir,'image7_R10_measurement_diagnostic_summary.csv'));
fprintf('Saved individual course files, results, and CSV summary in %s\n',outDir);
fprintf('A PLAUSIBLE numerical result is not proof of convergence.\n');
fprintf('Synthetic replacement trials are DIAGNOSTIC ONLY. Do not use their adjusted coordinates as final results.\n');

function H=fitPlanarH(x1,x2)
% Normalized DLT: x2 ~= H*x1; at least 4 paired planar targets.
assert(size(x1,1)>=4 && size(x1,1)==size(x2,1));
[xn,T1]=normPts(x1);[yn,T2]=normPts(x2);
n=size(xn,1);M=zeros(2*n,9);
for j=1:n
 x=xn(j,1);y=xn(j,2);u=yn(j,1);v=yn(j,2);
 M(2*j-1,:)=[-x -y -1 0 0 0 u*x u*y u];
 M(2*j,:)=[0 0 0 -x -y -1 v*x v*y v];
end
[~,~,VV]=svd(M,0);
H=T2\reshape(VV(:,end),3,3)'*T1;
end
function [p,T]=normPts(p)
 c=mean(p,1); s=sqrt(2)/max(mean(vecnorm(p-c,2,2)),eps);
 T=[s 0 -s*c(1);0 s -s*c(2);0 0 1];
 h=T*[p ones(size(p,1),1)]';p=h(1:2,:)';
end

%% ENGO 531 Lab 4 - SELECTABLE IMAGES: convergence-only trial
% Calls the course full_matlab_bundle.p. No data snooping.
% Uses Lab 1 mode 1: FIRST-image weighted EOP observations + one scale distance.
% The photogrammetric IOPs are fixed to Lab 3 estimates.
% Tie-point approximate XYZ are initialized from the photographed pod faces.
% One trial. All original geometry construction and P-code settings preserved.
function result = runOneTrial(phoFile,outDir,selectedImages,minImagesPerPoint,excludedTargets,trialTag,replaceImage,replaceXY)
result = struct('Trial',string(trialTag),'Excluded',strjoin(string(excludedTargets),','), ...
    'RMS_px',NaN,'MaxCameraShift_mm',NaN,'Finite',false,'Status',"ERROR",'Message',"");

%% SETTINGS
imageSize = [5712 4284];      % width,height in pixels
IOP = [2834.0 2156.7 3995.9]; % xp,yp,c (pixels); fixed
% Change this between runs. Examples: [1 2 3 4 5 6], [1:7], [1:9].
% Datum: the FIRST selected image gets six weighted EOP observations.
% Include Image 1 in comparisons if you want to preserve the same datum.
% selectedImages is supplied by the diagnostic driver
% minImagesPerPoint is supplied by the diagnostic driver       % Set to 2 for additional troubleshooting
% Change to true only to allow F1/F6 with two images when minimum is 3.
allowTwoImageScaleEndpoints = false;
scaleFrom = "F1"; scaleTo = "F6";
measuredDistance_mm = 2478;
sd_xy = 3.2;                  % image observations (pixels)
sd_dist = 5;                 % trial tape-measure SD (mm)
sd_firstEOP_XYZ = 1;         % first camera EOP weighted position SD (mm)
sd_firstEOP_wpk = 0.01;      % first camera angles SD (degrees)

% Original field-note EOPs with reversed phi: [image Xc Yc Zc phi kappa omega]
notes = [
 1 -2213  1711 -645  -90 0 90;
 2 -2279   561 -495  -75 0 90;
 3 -1619 -1540 -495  -45 0 90;
 4  -440 -2027 -495  -25 0 90;
 5  1239 -2545 -495    0 0 90;
 6  3486 -1652 -495   30 0 90;
 7  4368   -97 -495   60 0 90;
 8  4901   411 -495   90 0 90;
 9  4455  2686 -495  120 0 90];

%% INPUT FILES
imageTag = strjoin(string(selectedImages),'_');
prefix = fullfile(outDir, ['diag_' char(trialTag) '_images_' char(imageTag)]);

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
% Replace ONE R10 photo observation for controlled diagnosis. The predicted
% coordinate is based on other targets' homography, never the actual R10.
% Both R10 measurements stay in the adjustment; no rank-deficient orphan.
if replaceImage~=0
    q=find(pointName=="R10" & imageID==replaceImage);
    assert(numel(q)==1,'Expected exactly one R10 on selected Image %d',replaceImage);
    fprintf('DIAGNOSTIC ONLY: Image %d R10 from [%.3f %.3f] to [%.3f %.3f] px\n', ...
        replaceImage,xy(q,:),replaceXY);
    xy(q,:)=replaceXY;
end
% Exclude target IDs from ALL selected images, without editing the PHO file.
if ~isempty(excludedTargets)
    keepTarget = ~ismember(pointName,string(excludedTargets));
    pointName=pointName(keepTarget); imageID=imageID(keepTarget); xy=xy(keepTarget,:);
end

fprintf('Selected images: %s\n',mat2str(activeImages'));

%% RETAIN POINTS SEEN IN THE SELECTED SUBSET
[names0,~,g]=unique(pointName,'stable'); counts=accumarray(g,1);
if strcmp(trialTag,'baseline'), disp(table(names0,counts,'VariableNames',{'Target','ImageCount'})); end
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
[tf,scaleIdx]=ismember([scaleFrom;scaleTo],names);assert(all(tf));
D0=norm(XYZ0(scaleIdx(2),:)-XYZ0(scaleIdx(1),:));
writeRows([prefix '.dis'],[numericID(scaleIdx(1)),numericID(scaleIdx(2)),measuredDistance_mm], ...
    '%d %d %.6f\n');
% CFG: param_obs_mode=1 (first image EOP observations); iop_mode=0 (fixed)
cfg=[1;0;0.001;0.0001;0.01;sd_xy;sd_dist;sd_firstEOP_XYZ;sd_firstEOP_wpk];
writeRows([prefix '.cfg'],cfg,'%g\n');
if isfile([prefix '.con']),delete([prefix '.con']);end

nUnknown=6*numImages+3*numPoints;
nScalar=2*numPhoto+6+1;
fprintf('\nSELECTED-IMAGE NETWORK: %d images, %d points, %d photo measurements\n', ...
    numImages,numPoints,numPhoto);
fprintf('Unknowns %d, scalar observations %d, nominal DOF %d\n',nUnknown,nScalar,nScalar-nUnknown);
fprintf('Initial F1-F6 distance: %.3f mm; tape distance: %.3f mm\n',D0,measuredDistance_mm);
fprintf('Selected IDs: %s; first image %d defines EOP datum.\n', mat2str(activeImages'),activeImages(1));
fprintf('Field-note reversed phi; fixed calibrated IOP.\n');
fprintf('No automatic data snooping. Synthetic R10 replacements are diagnostics only. Files: %s.*\n',prefix);
if nScalar-nUnknown<=0
    warning('Nominal degrees of freedom are nonpositive; adjustment is unlikely to be identifiable.');
end

%% RUN COURSE BUNDLE ADJUSTMENT
assert(exist('full_matlab_bundle','file')==2 || exist('full_matlab_bundle','file')==6, ...
    'Add full_matlab_bundle.p to MATLAB path');
[eop,iop,con,tie,Cx,vr_img,vr_eop,vr_dis,vr_con]=full_matlab_bundle(prefix);
rmsXY=sqrt(mean(vr_img(:,1).^2));
result.RMS_px=rmsXY;
result.MaxCameraShift_mm=max(vecnorm(eop.XYZc-ext(:,3:5),2,2));
result.Finite=all(isfinite(eop.XYZc(:))) && all(isfinite(eop.wpk(:))) && ...
    all(isfinite(vr_img(:,1)));
% This is a gross-instability screen, NOT an independent convergence check.
if ~result.Finite || result.MaxCameraShift_mm>1e5 || rmsXY>100
    result.Status="DIVERGED/INVALID";
else
    result.Status="PLAUSIBLE--VERIFY MATLAB CONVERGENCE MESSAGE";
end

fprintf('Returned RMS = %.3f px (ONLY meaningful if MATLAB reports convergence)\n',rmsXY);
disp('Returned adjusted camera XYZ [mm]:'); disp(eop.XYZc);
disp('Returned adjusted camera omega/phi/kappa [degrees]:'); disp(eop.wpk);
if ~isempty(vr_dis),fprintf('Distance residual: %.3f mm\n',vr_dis(1,1));end
save([prefix '_results.mat'],'eop','iop','con','tie','Cx','vr_img', ...
    'vr_eop','vr_dis','vr_con','ext','IOP','cfg','names','numericID','XYZ0', ...
    'photoTable','tieTable','rmsXY','result');

end % runOneTrial

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
