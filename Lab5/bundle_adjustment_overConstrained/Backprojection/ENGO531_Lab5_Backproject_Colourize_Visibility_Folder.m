function ENGO531_Lab5_Backproject_Colourize_Visibility_Folder()
% ENGO 531 Lab 5 - scanner XYZ point cloud to RGB with face-aware visibility.
% MATLAB R2024b. Does NOT alter the bundle adjustment.
% Uses the verified collinearity convention from the earlier Lab 5 script.
%
% 1) Project validation against measured image points.
% 2) Fit left/front/right planes using the TLS picked point CSV.
% 3) Assign every TLS point to the closest suitable box-face plane.
% 4) Restrict candidate cameras by face; reject points behind nearer scanned
%    points using a coarse image-space Z-buffer; choose the most frontal view.
% 5) Save six-column XYZRGB for CloudCompare, and a diagnostics summary.
%
% CAUTION: Glass and occluded faces make visibility estimation imperfect.
% Inspect the output; vary face tolerance and depth tolerance as needed.
clc;
fprintf('\nENGO 531 LAB 5: FACE-AWARE TLS COLOURIZATION\n');

%% User-tunable settings
faceGroups.L = [1 2 3];
faceGroups.F = [4 5 6];
faceGroups.R = [7 8 9];
FACE_TOLERANCE_MM = 130;      % reject points far from the three fitted planes
DEPTH_TOLERANCE_MM = 35;     % allowable depth difference from nearest point
ZBUFFER_BIN_PX = 4;          % depth-buffer cell size in full-resolution pixels
MAX_REPROJECTION_RMS_PX = 5; % stop if camera convention is incompatible
NORMAL_MIN_COSINE = 0.08;    % reject near-tangent, back-facing camera views
MAX_INCIDENCE_WEIGHT = 0.65;% choose mostly frontal views, small centre penalty
PREVIEW_POINT_LIMIT = 45000;

%% Adjusted Lab 5 bundle
[fn,fd]=uigetfile('*.mat','Select converged Lab 5 SIX-control *_results.mat');
if isequal(fn,0),return;end
S=load(fullfile(fd,fn));
req={'eop','tie','con','tieRows','conRows','pho','ext'};
assert(all(isfield(S,req)), 'Bundle MAT missing needed fields.');
assert(all(isfield(S.eop,{'XYZc','wpk'})),'Bundle EOPs are incomplete.');
C=double(S.eop.XYZc);wpk=double(S.eop.wpk);
imageIDs=double(S.ext(:,1));nCam=numel(imageIDs);
assert(size(C,1)==nCam && size(wpk,1)==nCam,'EOP count mismatch');

intrFile=fullfile(fd,'lab5_TLS.int');
if ~isfile(intrFile)
    [inf,ind]=uigetfile('*.int','Select the .int used in Lab 5 bundle');
    if isequal(inf,0),return;end
    intrFile=fullfile(ind,inf);
end
intr=readmatrix(intrFile,'FileType','text');
IOP=double(intr(1,2:4));
assert(all(isfinite(IOP)) && IOP(3)>0,'Invalid IOP');
fprintf('Fixed IOP [xp yp c] = %.4f %.4f %.4f pixels\n',IOP);

%% Original photographs, automatically matched to image IDs
imageDir=uigetdir(fd,'Select folder containing Image1.jpg, Image2.jpg, ...');
if isequal(imageDir,0),return;end
images=cell(nCam,1);W=zeros(nCam,1);H=zeros(nCam,1);
imagePaths=strings(nCam,1);
fileList=dir(imageDir);
fileList=fileList(~[fileList.isdir]);
fileNames=string({fileList.name});
for k=1:nCam
    % Match exact Image<number> filename (case insensitive); accept JPG/PNG/TIF.
    % Leading zeros are also accepted (e.g., Image01.jpg for Image 1).
    matches=false(numel(fileNames),1);
    for fidx=1:numel(fileNames)
        tokens=regexpi(char(fileNames(fidx)), ...
            '^Image([0-9]+)\.(jpg|jpeg|png|tif|tiff)$','tokens','once');
        if ~isempty(tokens) && str2double(tokens{1})==imageIDs(k)
            matches(fidx)=true;
        end
    end
    candidates=fileNames(matches);
    assert(~isempty(candidates), ...
        'Missing Image%d photograph in %s. Expected Image%d.jpg (or PNG/TIF).', ...
        imageIDs(k),imageDir,imageIDs(k));
    assert(numel(candidates)==1, ...
        'Multiple photographs match Image%d: %s. Keep only one.', ...
        imageIDs(k),strjoin(candidates,', '));
    imagePaths(k)=string(fullfile(imageDir,candidates(1)));
    im=imread(char(imagePaths(k)));
    assert(size(im,3)>=3,'Image %d must contain RGB channels',imageIDs(k));
    if ~isa(im,'uint8'),im=im2uint8(im(:,:,1:3));end
    images{k}=im(:,:,1:3);H(k)=size(im,1);W(k)=size(im,2);
    fprintf('Image %d: %s (%d x %d)\n',imageIDs(k),candidates(1),W(k),H(k));
end

%% Validate projection geometry against image observations
pointIDs=[double(S.tieRows(:,1));double(S.conRows(:,1))];
pointXYZ=[double(S.tie.XYZ);double(S.con.XYZ)];
assert(numel(unique(pointIDs))==numel(pointIDs),'Duplicate adjusted point IDs');
pho=double(S.pho);
[found,ir]=ismember(pho(:,1),pointIDs);
assert(all(found),'Measured image point not found in adjusted XYZ');
obsXYZ=pointXYZ(ir,:);
depthSigns=zeros(nCam,1);allErr=[];
for k=1:nCam
    q=pho(:,2)==imageIDs(k);
    assert(any(q),'Missing PHO rows for Image %d',imageIDs(k));
    [pix,d]=projectXYZ(obsXYZ(q,:),C(k,:),wpk(k,:),IOP);
    err=hypot(pix(:,1)-pho(q,3),pix(:,2)-pho(q,4));
    allErr=[allErr;err]; %#ok<AGROW>
    depthSigns(k)=sign(median(d,'omitnan'));
    assert(depthSigns(k)~=0,'Cannot establish camera depth sign');
    fprintf('Image %d: radial RMS %.3f px\n',imageIDs(k),sqrt(mean(err.^2,'omitnan')));
end
rmsPixel=sqrt(mean(allErr.^2,'omitnan'));
fprintf('Overall radial reprojection RMS = %.3f px\n',rmsPixel);
assert(rmsPixel<=MAX_REPROJECTION_RMS_PX, ...
    'Reprojection error is too large: check EOP/IOP or images');

%% TLS picked target coordinates for face plane fitting
[tfn,tfd]=uigetfile('*.csv','Choose target_coordinates CSV (point_id,X,Y,Z)');
if isequal(tfn,0),return;end
T=readtable(fullfile(tfd,tfn),'TextType','string');
assert(all(ismember({'point_id','X','Y','Z'},T.Properties.VariableNames)), ...
    'Expected TLS target CSV: point_id,X,Y,Z,cloud_index');
pointName=upper(strtrim(string(T.point_id)));
tlsTargetXYZ=1000*double([T.X T.Y T.Z]); % CSV known metres
assert(all(isfinite(tlsTargetXYZ(:))),'Invalid TLS target XYZ');
faceLetters='LFR';
planeCentres=zeros(3,3);planeNormals=zeros(3,3);
planeRMS=zeros(3,1);
for j=1:3
    rows=startsWith(pointName,string(faceLetters(j)));
    P=tlsTargetXYZ(rows,:);
    assert(size(P,1)>=3,'Need at least 3 picked TLS targets for face %s',faceLetters(j));
    centre=mean(P,1);
    [~,~,V]=svd(P-centre,0);
    normal=V(:,end)';normal=normal/norm(normal);
    % Orient the face normal toward the cameras nominally viewing that face.
    camIDs=faceGroups.(faceLetters(j));
    mask=ismember(imageIDs,camIDs);
    assert(any(mask),'No cameras assigned to face %s',faceLetters(j));
    camMean=mean(C(mask,:),1);
    if dot(normal,camMean-centre)<0,normal=-normal;end
    planeCentres(j,:)=centre;planeNormals(j,:)=normal;
    dist=(P-centre)*normal';
    planeRMS(j)=sqrt(mean(dist.^2));
    fprintf('%s plane: %d targets; fit RMS %.2f mm; normal [%.3f %.3f %.3f]\n', ...
        faceLetters(j),size(P,1),planeRMS(j),normal);
end
if any(planeRMS>25)
    warning('One or more face planes fit poorly (>25 mm). Check the TLS target face IDs.');
end

%% Cropped point cloud in the scanner datum
[pcf,pcd]=uigetfile({'*.txt;*.xyz;*.asc;*.csv','Point cloud'}, ...
    'Select CROPPED raw TLS cloud (XYZ or XYZRGB)');
if isequal(pcf,0),return;end
cloudPath=fullfile(pcd,pcf);
unitChoice=menu('Raw TLS point cloud XYZ units?','Metres','Millimetres');
if unitChoice==0,return;end
if unitChoice==1,mmPerCloudUnit=1000;else,mmPerCloudUnit=1;end
[XYZraw,nSkipped]=loadXYZ(cloudPath);
N=size(XYZraw,1);
assert(N>0,'No valid XYZ points in input');
XYZ=XYZraw*mmPerCloudUnit;
fprintf('Loaded %d raw TLS points; skipped %d header/invalid lines.\n',N,nSkipped);

%% Classify point to its nearest fitted box plane
planeDistances=zeros(N,3,'single');
for j=1:3
    planeDistances(:,j)=single(abs((XYZ-planeCentres(j,:))*planeNormals(j,:)'));
end
[distanceToFace,faceNum]=min(planeDistances,[],2);
validFace=distanceToFace<=FACE_TOLERANCE_MM;
fprintf('Face classification: L=%d F=%d R=%d; outside tolerance=%d\n', ...
    nnz(validFace & faceNum==1),nnz(validFace & faceNum==2), ...
    nnz(validFace & faceNum==3),nnz(~validFace));

%% Per-camera Z-buffer computed from ALL points (not just selected face)
% This rejects points that project behind another observed TLS surface.
% A coarse grid reduces memory. It is an approximation, especially for glass.
Z=cell(nCam,1);BW=zeros(nCam,1);BH=zeros(nCam,1);
for k=1:nCam
    [pix,depth]=projectXYZ(XYZ,C(k,:),wpk(k,:),IOP);
    dep=depth*depthSigns(k);
    [bx,by,ok,bw,bh]=bufferCoordinates(pix,dep,W(k),H(k),ZBUFFER_BIN_PX);
    BW(k)=bw;BH(k)=bh;
    ids=sub2ind([bh,bw],by(ok),bx(ok));
    closeDepth = accumarray(ids,dep(ok),[bw*bh 1], ...
    @min,realmax('double'));
    Z{k}=single(closeDepth);
    fprintf('Depth buffer for Image %d complete (%d points in frame).\n', ...
        imageIDs(k),nnz(ok));
end

%% Choose best image by face, visibility, frontality and image centrality
bestScore = repmat(realmax('double'), N, 1);
bestCamera=zeros(N,1,'uint8');
RGB=zeros(N,3,'uint8');
acceptedByImage=zeros(nCam,1);
for k=1:nCam
    % Restrict candidate images to expected face groups.
    eligible=false(N,1);
    for j=1:3
        if ismember(imageIDs(k),faceGroups.(faceLetters(j)))
            eligible=eligible | (faceNum==j & validFace);
        end
    end
    if ~any(eligible),continue;end
    indices=find(eligible);
    P=XYZ(indices,:);
    [pix,depth]=projectXYZ(P,C(k,:),wpk(k,:),IOP);
    dep=depth*depthSigns(k);
    [bx,by,inFrame,~,~]=bufferCoordinates(pix,dep,W(k),H(k),ZBUFFER_BIN_PX);
    if ~any(inFrame),continue;end
    indLocal=find(inFrame);
    pix=pix(indLocal,:);P=P(indLocal,:);dep=dep(indLocal);
    globalIdx=indices(indLocal);
    id=sub2ind([BH(k),BW(k)],by(inFrame),bx(inFrame));
    visible=dep<=double(Z{k}(id))+DEPTH_TOLERANCE_MM;
    if ~any(visible),continue;end
    pix=pix(visible,:);P=P(visible,:);globalIdx=globalIdx(visible);
    % Frontal incidence: angle between surface normal and camera-point ray.
    ray=C(k,:)-P;
    rayNorm=sqrt(sum(ray.^2,2));
    norms=planeNormals(faceNum(globalIdx),:);
    cosine=sum(ray.*norms,2)./max(rayNorm,eps);
    front=cosine>=NORMAL_MIN_COSINE;
    if ~any(front),continue;end
    pix=pix(front,:);globalIdx=globalIdx(front);cosine=cosine(front);
    % Penalize strongly oblique viewing and image periphery.
    normalizedRadius=((pix(:,1)-(W(k)+1)/2)/(W(k)/2)).^2 + ...
                     ((pix(:,2)-(H(k)+1)/2)/(H(k)/2)).^2;
    score=MAX_INCIDENCE_WEIGHT*(1-cosine).^2+(1-MAX_INCIDENCE_WEIGHT)*normalizedRadius;
    take=score<bestScore(globalIdx);
    if ~any(take),continue;end
    ids=globalIdx(take);
    pp=pix(take,:);
    col=round(pp(:,1));row=round(pp(:,2));
    linear=sub2ind([H(k),W(k)],row,col);
    im=images{k};
    for ch=1:3
        layer=im(:,:,ch);
        RGB(ids,ch)=layer(linear);
    end
    bestScore(ids)=score(take);
    bestCamera(ids)=uint8(k);
    fprintf('Image %d: currently assigned %d points.\n',imageIDs(k),nnz(bestCamera==k));
end
keep=bestCamera>0;

%% Export X Y Z R G B, without a header, retaining cloud units
[outfn,outfd]=uiputfile('*.txt','Save FACE-AWARE XYZRGB', ...
    fullfile(pcd,'lab5_TLS_colourized_VISIBILITY_XYZRGB.txt'));
if isequal(outfn,0),return;end
outPath=fullfile(outfd,outfn);
fid=fopen(outPath,'wt');assert(fid>0,'Unable to create XYZRGB output');
fcloseGuard=onCleanup(@()fclose(fid));
block=25000;validIdx=find(keep);
for a=1:block:numel(validIdx)
    sel=validIdx(a:min(a+block-1,numel(validIdx)));
    values=[XYZraw(sel,:),double(RGB(sel,:))];
    fprintf(fid,'%.8f %.8f %.8f %d %d %d\n',values');
end
clear fcloseGuard
logPath=fullfile(outfd,'lab5_visibility_colourization_report.txt');
logFID=fopen(logPath,'wt');
if logFID>0
    fprintf(logFID,'ENGO 531 Lab 5 face-aware TLS colourization\n');
    fprintf(logFID,'Bundle: %s\nPoint cloud: %s\nTarget CSV: %s\n',fullfile(fd,fn),cloudPath,fullfile(tfd,tfn));
    fprintf(logFID,'Reprojection radial RMS %.4f px\n',rmsPixel);
    fprintf(logFID,'Plane fit RMS (L F R) mm: %.4f %.4f %.4f\n',planeRMS);
    fprintf(logFID,'Face tolerance %.1f mm; Zbuffer size %d px; depth tolerance %.1f mm\n', ...
        FACE_TOLERANCE_MM,ZBUFFER_BIN_PX,DEPTH_TOLERANCE_MM);
    fprintf(logFID,'Input points %d; saved points %d; excluded points %d\n',N,nnz(keep),nnz(~keep));
    for k=1:nCam
        count=nnz(bestCamera==k);
        fprintf(logFID,'Image %d: %d assigned points; %s\n',imageIDs(k),count,imagePaths(k));
    end
    fclose(logFID);
end
fprintf('\nSAVED: %d/%d (%.1f%%) TLS points colourized.\n', ...
    nnz(keep),N,100*nnz(keep)/N);
fprintf('File: %s\nLog: %s\n',outPath,logPath);

% Optional MATLAB preview for quick inspection.
if nnz(keep)>0
    idx=find(keep);
    idx=idx(round(linspace(1,numel(idx),min(numel(idx),PREVIEW_POINT_LIMIT))));
    figure('Name','Lab 5 - visible point RGB preview');
    scatter3(XYZraw(idx,1),XYZraw(idx,2),XYZraw(idx,3),2,double(RGB(idx,:))/255,'filled');
    axis equal;grid on;view(3);
    title(sprintf('Face-aware RGB preview: %.1f%% coloured',100*nnz(keep)/N));
    xlabel('Scanner X');ylabel('Scanner Y');zlabel('Scanner Z');
end
fprintf('Inspect in CloudCompare. Remaining seams may be caused by glass, TLS pick error or camera offsets.\n');
end

function [pixel,depth]=projectXYZ(XYZ,cameraXYZ,wpk,IOP)
    a=deg2rad(wpk(1));b=deg2rad(wpk(2));g=deg2rad(wpk(3));
    Rx=[1 0 0;0 cos(a) -sin(a);0 sin(a) cos(a)];
    Ry=[cos(b) 0 sin(b);0 1 0;-sin(b) 0 cos(b)];
    Rz=[cos(g) -sin(g) 0;sin(g) cos(g) 0;0 0 1];
    M=(Rx*Ry*Rz)';
    q=(M*(XYZ-cameraXYZ)')';
    depth=q(:,3);
    d=depth;d(abs(d)<1e-10)=NaN;
    pixel=[IOP(1)-IOP(3)*q(:,1)./d,IOP(2)+IOP(3)*q(:,2)./d];
end

function [bx,by,valid,bw,bh]=bufferCoordinates(pix,dep,width,height,bin)
    bw=ceil(width/bin);bh=ceil(height/bin);
    valid=isfinite(pix(:,1)) & isfinite(pix(:,2)) & isfinite(dep) & ...
        dep>0 & pix(:,1)>=1 & pix(:,1)<=width & ...
        pix(:,2)>=1 & pix(:,2)<=height;
    bx=min(bw,max(1,floor((pix(:,1)-1)/bin)+1));
    by=min(bh,max(1,floor((pix(:,2)-1)/bin)+1));
end

function [XYZ,skipped]=loadXYZ(path)
    fid=fopen(path,'rt');assert(fid>0,'Cannot open cloud');
    cleaner=onCleanup(@()fclose(fid));
    chunk=50000;parts=cell(0,1);skipped=0;
    while true
        A=zeros(chunk,3);n=0;endOfFile=false;
        for q=1:chunk
            ln=fgetl(fid);
            if ~ischar(ln),endOfFile=true;break;end
            ln=strtrim(ln);
            if isempty(ln) || startsWith(ln,'#') || startsWith(ln,'//')
                skipped=skipped+1;continue;
            end
            vals=sscanf(regexprep(ln,'[,;]',' '),'%f');
            if numel(vals)>=3 && all(isfinite(vals(1:3)))
                n=n+1;A(n,:)=vals(1:3)';
            else
                skipped=skipped+1;
            end
        end
        if n>0,parts{end+1,1}=A(1:n,:);end %#ok<AGROW>
        if endOfFile,break;end
    end
    if isempty(parts),XYZ=zeros(0,3);else,XYZ=vertcat(parts{:});end
end
