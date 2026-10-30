%% ENGO 531 Lab 4: RIGHT-WALL-ONLY bundle adjustment (Images 7, 8, 9)
% Standalone course full_matlab_bundle.p run; no other lab scripts required.
% Only R targets are included. Six weighted EOP observations for Image 7
% define the datum. The measured R4--R1 distance (2701 mm) sets scale.
% IOP are FIXED to the calibrated Lab 3 camera parameters.
% All XYZ quantities are in mm; image coordinates and IOPs are pixels.
clear; clc; close all;

%% SETTINGS
selectedImages = [7 8 9];
minImagesPerPoint = 2;             % all three-image R targets are preferred
IOP = [2834.0 2156.7 3995.9];     % xp, yp, c [pixels] calibrated
rightWallDepth_mm = 2701;         % R4 to R1, physically measured
rightWallHeight_mm = 2213;        % approximate, for XYZ initialization ONLY
rightWallX_mm = 2543;             % approximate X, not control
scaleFrom = "R4"; scaleTo = "R1";
measuredDistance_mm = 2701;       % R4--R1 measured 2.701 m
sd_xy = 3.2;                      % image-coordinate SD [px]
sd_dist = 5;                      % ASSUMED tape-measure SD [mm]; edit if known
sd_firstEOP_XYZ = 1;              % first camera weighted datum SD [mm]
sd_firstEOP_wpk = 0.01;           % first camera datum SD [degrees]

% Field notes: [image Xc Yc Zc phi kappa omega] (in the same frame as prior scripts)
notes = [
    7 4368  -97 -495   60 0 90;
    8 4901  411 -495   90 0 90;
    9 4455 2686 -495  120 0 90];

%% SELECT SOURCE PHO AND OUTPUT FOLDER
[fn,fd] = uigetfile({'*.PHO;*.pho','PHO point measurements'}, ...
    'Select calibration(4).PHO (or latest complete PHO)');
if isequal(fn,0), return; end
src = fullfile(fd,fn);
outDir = uigetdir(fd,'Select output folder for right-wall 7-8-9 test');
if isequal(outDir,0), return; end
prefix = fullfile(outDir,'engo531_right_wall_images_7_8_9');

fid = fopen(src,'r'); assert(fid>=0,'Cannot open PHO');
C = textscan(fid,'%s %f %f %f'); fclose(fid);
allNames = string(C{1}); allImages = C{2}; allXY = [C{3} C{4}];
assert(~isempty(allNames),'PHO is empty.');
assert(all(isfinite([allImages;allXY(:)])),'PHO has invalid numbers.');
assert(all(ismember(selectedImages,unique(allImages))), 'PHO missing selected images.');

% Select ONLY right-wall targets from Images 7, 8 and 9.
sel = startsWith(allNames,'R') & ismember(allImages,selectedImages);
pName = allNames(sel); imID = allImages(sel); xy = allXY(sel,:);
assert(~isempty(pName),'No right-wall observations found.');
assert(numel(unique(pName+"_"+string(imID))) == numel(pName), ...
    'A target was measured twice in the same image.');
[names0,~,g] = unique(pName,'stable');
counts0 = accumarray(g,1);
disp(table(names0,counts0,'VariableNames',{'Target','ImageCount'}));
keep = names0(counts0 >= minImagesPerPoint);
sel2 = ismember(pName,keep);
pName = pName(sel2); imID = imID(sel2); xy = xy(sel2,:);
[names,~,pointIdx] = unique(pName,'stable');
assert(all(ismember([scaleFrom;scaleTo],names)), ...
    'R4 and R1 must both be present to use their measured distance.');
assert(all(ismember(selectedImages,unique(imID))), ...
    'Some selected image has no retained observations.');

% Check connectivity by shared targets.
seen = false(numel(selectedImages),1); seen(1) = true;
changed = true;
while changed
    changed = false;
    for a = 1:numel(selectedImages)
        if ~seen(a), continue; end
        A = unique(pName(imID == selectedImages(a)));
        for b = 1:numel(selectedImages)
            if ~seen(b) && any(ismember(A,pName(imID == selectedImages(b))))
                seen(b) = true; changed = true;
            end
        end
    end
end
assert(all(seen),'Right-wall image network is disconnected.');

%% INITIAL XYZ FROM THE RIGHT-WALL FOUR CORNERS IN IMAGE 8
% Corner mapping: R4=(Y=0,Z=0), R1=(Y=2701,Z=0),
% R8=(Y=2701,Z=-2213), R5=(Y=0,Z=-2213).
cornerNames = ["R4" "R1" "R8" "R5"];
cornerUV = [0 0; rightWallDepth_mm 0; ...
    rightWallDepth_mm -rightWallHeight_mm; 0 -rightWallHeight_mm];
H = localFaceH(cornerNames,8,allNames,allImages,allXY,cornerUV);

nPts = numel(names); XYZ0 = zeros(nPts,3);
for j=1:nPts
    q = find(allNames==names(j) & allImages==8,1);
    assert(~isempty(q),'Point %s is absent from Image 8; cannot initialize.',names(j));
    uv = localApplyH(H,allXY(q,:));
    XYZ0(j,:) = [rightWallX_mm uv(1) uv(2)];
end
id = (1001:(1000+nPts))';
idMap = table(names,id,'VariableNames',{'OriginalPointID','BundlePointID'});
writetable(idMap,[prefix '_target_ID_map.csv']);
writetable(table(names,XYZ0(:,1),XYZ0(:,2),XYZ0(:,3), ...
    'VariableNames',{'Target','X_mm','Y_mm','Z_mm'}), ...
    [prefix '_initial_XYZ.csv']);

%% BUILD AND RUN THE COURSE BUNDLE
% .ext format: image camera X Y Z omega phi kappa
ext = [notes(:,1),99*ones(size(notes,1),1),notes(:,2:4), ...
       notes(:,7),notes(:,5),notes(:,6)];
phoTable = [id(pointIdx),imID,xy];
tieTable = [id,XYZ0];
localWrite([prefix '.pho'],phoTable,'%d %d %.4f %.4f\n');
localWrite([prefix '.tie'],tieTable,'%d %.6f %.6f %.6f\n');
localWrite([prefix '.ext'],ext,'%d %d %.6f %.6f %.6f %.6f %.6f %.6f\n');
localWrite([prefix '.int'],[99 IOP -1],'%d %.6f %.6f %.6f %d\n');
[~,jScale] = ismember([scaleFrom;scaleTo],names);
localWrite([prefix '.dis'],[id(jScale(1)),id(jScale(2)),measuredDistance_mm], ...
    '%d %d %.6f\n');
% .cfg: param_obs_mode=1, iop_mode=0, convergence tolerances, observation SDs
cfg = [1;0;0.001;0.0001;0.01;sd_xy;sd_dist;sd_firstEOP_XYZ;sd_firstEOP_wpk];
localWrite([prefix '.cfg'],cfg,'%g\n');
if isfile([prefix '.con']), delete([prefix '.con']); end

nImages = numel(selectedImages); nPhoto = size(phoTable,1);
nUnknown = 6*nImages + 3*nPts; nObserv = 2*nPhoto + 6 + 1;
fprintf('\nRIGHT WALL: images %s, %d points, %d photo measurements\n', ...
    mat2str(selectedImages),nPts,nPhoto);
fprintf('Unknowns = %d, scalar observations = %d, nominal DOF = %d\n', ...
    nUnknown,nObserv,nObserv-nUnknown);
fprintf('Initial R4--R1 = %.3f mm; measured = %.3f mm\n', ...
    norm(XYZ0(jScale(1),:)-XYZ0(jScale(2),:)),measuredDistance_mm);
fprintf('Datum: weighted EOPs of Image 7. IOPs fixed (calibrated).\n');
fprintf('Files saved as %s.*\n',prefix);
assert(nObserv>nUnknown,'Not enough observations for an overdetermined network.');
assert(exist('full_matlab_bundle','file')==2 || ...
       exist('full_matlab_bundle','file')==6, ...
       'Add the course full_matlab_bundle.p to your MATLAB path.');

[eop,iop,con,tie,Cx,vr_img,vr_eop,vr_dis,vr_con] = full_matlab_bundle(prefix);
rmsXY = sqrt(mean(vr_img(:,1).^2));
fprintf('\nReturned RMS: %.3f px (VALID ONLY IF COURSE REPORTS CONVERGENCE)\n',rmsXY);
disp('Adjusted camera XYZ [mm]:'); disp(eop.XYZc);
disp('Adjusted omega/phi/kappa [degrees]:'); disp(eop.wpk);
if ~isempty(vr_dis), fprintf('R4--R1 distance residual = %.4f mm\n',vr_dis(1,1)); end
save([prefix '_results.mat'], 'eop','iop','con','tie','Cx', ...
    'vr_img','vr_eop','vr_dis','vr_con','names','id','idMap', ...
    'XYZ0','ext','cfg','IOP','phoTable','tieTable','rmsXY');

% Plot only if the user confirms the course routine converged.
answ = input('Did MATLAB explicitly report Adjustment converged? [y/N]: ','s');
if strcmpi(strtrim(answ),'y')
    try
        P = tie.XYZ;
        figure('Color','w');
        scatter3(P(:,1),P(:,2),P(:,3),60,'filled'); hold on;
        text(P(:,1),P(:,2),P(:,3),"  "+names);
        scatter3(eop.XYZc(:,1),eop.XYZc(:,2),eop.XYZc(:,3), ...
            90,'^','filled');
        for k=1:nImages
            text(eop.XYZc(k,1),eop.XYZc(k,2),eop.XYZc(k,3), ...
                sprintf('  Image %d',selectedImages(k)));
        end
        axis equal; grid on; rotate3d on;
        xlabel('X [mm]');ylabel('Y [mm]');zlabel('Z [mm]');
        title('Right-wall bundle: adjusted tie points and cameras');
        legend('R tie points','Camera stations','Location','best');
        exportgraphics(gcf,[prefix '_adjusted_3D.png'],'Resolution',200);
    catch ME
        warning('Adjustment returned, but plot failed: %s',ME.message);
    end
else
    fprintf('Plot suppressed because convergence was not confirmed.\n');
end

%% LOCAL HELPERS
function localWrite(path,A,fmt)
    f=fopen(path,'w'); assert(f>=0,'Cannot write %s',path);
    clean=onCleanup(@()fclose(f));
    for r=1:size(A,1),fprintf(f,fmt,A(r,:));end
end
function H=localFaceH(corners,img,allNames,allImgs,allXY,uv)
    A=zeros(8,8); b=zeros(8,1);
    for t=1:4
        ix=find(allNames==corners(t) & allImgs==img,1);
        assert(~isempty(ix),'Corner %s not found in Image %d',corners(t),img);
        x=allXY(ix,1);y=allXY(ix,2);u=uv(t,1);v=uv(t,2);
        A(2*t-1,:)=[x y 1 0 0 0 -u*x -u*y];
        A(2*t,:)=[0 0 0 x y 1 -v*x -v*y];
        b(2*t-1:2*t)=[u;v];
    end
    assert(rcond(A)>1e-14,'Ill-conditioned right-wall homography.');
    h=A\b;
    H=[h(1) h(2) h(3);h(4) h(5) h(6);h(7) h(8) 1];
end
function uv=localApplyH(H,xy)
    q=H*[xy(:);1];
    assert(abs(q(3))>1e-12,'Invalid homography transformation.');
    uv=(q(1:2)/q(3))';
end
