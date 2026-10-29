function ENGO531_Lab4_resection_warmstart()
% ENGO 531 - Six-image warm start + Image 7 planar resection diagnostic.
% Requires course full_matlab_bundle.p in MATLAB path.
% DOES NOT assume that R targets have 3D coordinates from a 1:6 adjustment.
% The R XYZ used in resection come from the approximate right-face mapping.
clc;
[fn,fd]=uigetfile({'*.PHO;*.pho','Photograph observations (*.PHO)'}, ...
    'Select calibration(4).PHO');
if isequal(fn,0),return;end
phoFile=fullfile(fd,fn);
outDir=uigetdir(fd,'Select output folder');
if isequal(outDir,0),return;end
fprintf('\n================ SIX-IMAGE REFERENCE ================\n');
base=runSubset(1:6,phoFile,outDir,[]);
fprintf('\nIMPORTANT: Confirm above says Adjustment converged BEFORE trusting warm start.\n');
reply=input('Did MATLAB report convergence for Images 1-6? Enter y to proceed: ','s');
if ~strcmpi(strtrim(reply),'y')
    fprintf('Stopped; no seven-image trial run.\n');return;
end
fprintf('\n================ SEVEN-IMAGE WARM START ================\n');
trial=runSubset(1:7,phoFile,outDir,base); %#ok<NASGU>
fprintf('\nReview image7_resection_residuals.csv and the bundle convergence message.\n');
end

function result=runSubset(imageSelection,phoFile,outDir,warm)
%% SETTINGS
imageSize = [5712 4284];      % width,height in pixels
IOP = [2834.0 2156.7 3995.9]; % xp,yp,c (pixels); fixed
% Change this between runs. Examples: [1 2 3 4 5 6], [1:7], [1:9].
% Datum: the FIRST selected image gets six weighted EOP observations.
% Include Image 1 in comparisons if you want to preserve the same datum.
selectedImages = imageSelection;
minImagesPerPoint = 2;       % Set to 2 for additional troubleshooting
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
prefix = fullfile(outDir, ['warm_images_' char(imageTag)]);

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

%% INPUT EXT: [image, camera, Xc Yc Zc omega phi kappa]
[found,rows]=ismember(activeImages,notes(:,1));assert(all(found));
chosenNotes=notes(rows,:);
ext=[chosenNotes(:,1),99*ones(numImages,1),chosenNotes(:,2:4), ...
     chosenNotes(:,7),chosenNotes(:,5),chosenNotes(:,6)];

%% APPLY SIX-IMAGE RESULTS AS INITIAL VALUES, THEN RESECT IMAGE 7
if ~isempty(warm)
    assert(isequal(selectedImages,1:7),'Warm start is intended for Images 1:7');
    assert(isfield(warm.tie,'XYZ') && size(warm.tie.XYZ,2)==3, ...
        'Six-image result has no usable adjusted tie XYZ');
    for jj=1:numPoints
        old=find(warm.names==names(jj),1);
        if ~isempty(old)
            XYZ0(jj,:)=warm.tie.XYZ(old,:);
        end
    end
    for ii=1:6
        ext(ii,3:5)=warm.eop.XYZc(ii,:);
        ext(ii,6:8)=warm.eop.wpk(ii,:);
    end
    % Right-wall XYZ remain approximate: R targets only appear in Image 6
    % during the six-image adjustment. Resection is NOT from measured 3D controls.
    ii7=find(imageID==7);
    rightSel=startsWith(pointName(ii7),'R');
    ii7=ii7(rightSel);
    [present,pt]=ismember(pointName(ii7),names);
    ii7=ii7(present);pt=pt(present);
    assert(numel(ii7)>=6,'Need at least six Image 7 right-wall observations');
    obs7=xy(ii7,:); object7=XYZ0(pt,:); labels7=pointName(ii7);
    guess=ext(ext(:,1)==7,3:8);
    [bestPose,bestErr,bestMode]=resectImage7(object7,obs7,guess,IOP);
    ext(ext(:,1)==7,3:8)=bestPose;
    fprintf('\nImage 7 resection projection convention: %s\n',bestMode);
    fprintf('Initial EOP7 XYZ/wpk: '); fprintf('%.3f ',guess);fprintf('\n');
    fprintf('Resected EOP7 XYZ/wpk: ');fprintf('%.3f ',bestPose);fprintf('\n');
    prefit=table(labels7,bestErr(:,1),bestErr(:,2), ...
        hypot(bestErr(:,1),bestErr(:,2)), ...
        'VariableNames',{'Target','ResidualX_px','ResidualY_px','RadialResidual_px'});
    prefit=sortrows(prefit,'RadialResidual_px','descend');
    disp(prefit);
    writetable(prefit,fullfile(outDir,'image7_resection_residuals.csv'));
    fprintf('Image 7 resection RMS radial: %.2f px, median %.2f px\n', ...
        sqrt(mean(prefit.RadialResidual_px.^2)),median(prefit.RadialResidual_px));
end

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
fprintf('No data snooping. Files: %s.*\n',prefix);
if nScalar-nUnknown<=0
    warning('Nominal degrees of freedom are nonpositive; adjustment is unlikely to be identifiable.');
end

%% RUN COURSE BUNDLE ADJUSTMENT
assert(exist('full_matlab_bundle','file')==2 || exist('full_matlab_bundle','file')==6, ...
    'Add full_matlab_bundle.p to MATLAB path');
[eop,iop,con,tie,Cx,vr_img,vr_eop,vr_dis,vr_con]=full_matlab_bundle(prefix);
rmsXY=sqrt(mean(vr_img(:,1).^2));
fprintf('Returned RMS = %.3f px (ONLY meaningful if MATLAB reports convergence)\n',rmsXY);
disp('Returned adjusted camera XYZ [mm]:'); disp(eop.XYZc);
disp('Returned adjusted camera omega/phi/kappa [degrees]:'); disp(eop.wpk);
if ~isempty(vr_dis),fprintf('Distance residual: %.3f mm\n',vr_dis(1,1));end
save([prefix '_results.mat'],'eop','iop','con','tie','Cx','vr_img', ...
    'vr_eop','vr_dis','vr_con','ext','IOP','cfg','names','numericID','XYZ0', ...
    'photoTable','tieTable','rmsXY','pointName','imageID','xy');


result=struct('eop',eop,'tie',tie,'names',names,'prefix',prefix, ...
    'rmsXY',rmsXY,'ext',ext,'imageID',imageID,'pointName',pointName);

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

function [poseOut,resOut,modeOut]=resectImage7(X,obs,guess,IOP)
% Freely initialized robust 6-parameter pose; two Y-image sign conventions
% tested because the course implementation is proprietary.
% Orientation model MUST be checked against course collinearity convention.
% Coordinates in km-scaled units for stable optimization.
starts=[guess; guess(1:3),guess(4:6)+[0 0 180]; ...
        guess(1:3),guess(4:6)+[0 180 0]; ...
        guess(1:3),[67.23 44.50 26.37]; ...
        guess(1:3),[90 60 0]];
opts=optimset('Display','off','MaxIter',5000,'MaxFunEvals',12000, ...
    'TolX',1e-8,'TolFun',1e-8);
best=Inf; poseOut=guess;resOut=NaN(size(obs));modeOut='';
for ySign=[1,-1]
    for j=1:size(starts,1)
        s=starts(j,:);
        q0=[s(1:3)/1000,s(4:6)/60];
        objective=@(q) reprojectionObjective(q,X,obs,IOP,ySign);
        [q,f]=fminsearch(objective,q0,opts);
        if f<best
            best=f;poseOut=[q(1:3)*1000,q(4:6)*60];
            resOut=projectResidual(poseOut,X,obs,IOP,ySign);
            modeOut=sprintf('Y image sign %d; Rz(-k)*Ry(phi)*Rx(-omega)',ySign);
        end
    end
end
end
function f=reprojectionObjective(q,X,obs,IOP,ySign)
pose=[q(1:3)*1000,q(4:6)*60];
r=projectResidual(pose,X,obs,IOP,ySign);
d=hypot(r(:,1),r(:,2));
% Smooth robust loss (Cauchy) limits leverage of suspected mislabeled target
f=sum(25^2*log1p((d/25).^2));
if ~isfinite(f),f=1e20;end
end
function residual=projectResidual(pose,X,obs,IOP,ySign)
w=deg2rad(pose(4));p=deg2rad(pose(5));k=deg2rad(pose(6));
Rx=[1 0 0;0 cos(w) -sin(w);0 sin(w) cos(w)];
Ry=[cos(p) 0 sin(p);0 1 0;-sin(p) 0 cos(p)];
Rz=[cos(k) -sin(k) 0;sin(k) cos(k) 0;0 0 1];
M=Rz'*Ry*Rx';
P=M*(X-pose(1:3))';
if any(abs(P(3,:))<1e-4)
    residual=ones(size(obs))*1e8;return;
end
pred=[IOP(1)-IOP(3)*P(1,:)'./P(3,:)', ...
      IOP(2)+ySign*IOP(3)*P(2,:)'./P(3,:)'];
residual=pred-obs;
if any(~isfinite(residual(:))),residual=ones(size(obs))*1e8;end
end
