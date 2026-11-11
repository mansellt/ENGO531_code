function ENGO531_Lab5_Backproject_Colourize_TLS()
% ENGO 531 Lab 5 - Colourize a cropped TLS XYZ cloud using adjusted camera EOPs.
%
% Inputs:
%   - CONVERGED six-control lab5_TLS_results.mat from Lab 5 Mode 0
%   - lab5_TLS.int (in same folder as the .mat, with fixed IOP [xp yp c])
%   - cropped CloudCompare text export: X Y Z [optional other columns]
%   - the nine ORIGINAL camera photographs (one selected for each image ID)
%
% Output: six-column SPACE-separated X Y Z R G B text, read in CloudCompare.
% Uses nearest-neighbour colour, choosing the image with the most central
% valid projection. This is a simple documented selection method, NOT a
% visibility/occlusion solution.
%
% The EOP rotation and pixel convention match the geometry-tested Lab 4
% MATLAB code: M = (Rx(omega)*Ry(phi)*Rz(kappa))';
% x = xp - c*Qc_x/Qc_z, y = yp + c*Qc_y/Qc_z.
% Coordinate inputs to projection must be in MILLIMETRES.
%
% IMPORTANT: This file is written for MATLAB R2024b or later.
clc;
fprintf('\nENGO 531 LAB 5 - TLS RGB BACKPROJECTION\n');

%% Choose adjustment
[fn,fd] = uigetfile('*.mat','Choose CONVERGED six-control Lab 5 results MAT');
if isequal(fn,0), return; end
S = load(fullfile(fd,fn));
required = {'eop','tie','con','tieRows','conRows','pho','ext'};
assert(all(isfield(S,required)), 'The MAT file lacks required fields: %s',strjoin(required,', '));
assert(isfield(S.eop,'XYZc') && isfield(S.eop,'wpk'), 'Missing adjusted EOPs');

% Read fixed IOPs from EXACT input used in the Lab 5 adjustment.
intFile = fullfile(fd,'lab5_TLS.int');
if ~isfile(intFile)
    [ifn,ifd] = uigetfile('*.int','Select .int actually used in this Lab 5 adjustment');
    if isequal(ifn,0),return;end
    intFile = fullfile(ifd,ifn);
end
intr = readmatrix(intFile,'FileType','text');
assert(size(intr,2)>=4 && size(intr,1)>=1, 'Expected .int: cameraID xp yp c ...');
camID = intr(1,1);
IOP = double(intr(1,2:4));  % [xp,yp,c], pixels
assert(all(isfinite(IOP)) && IOP(3)>0,'Invalid interior orientation');
fprintf('Using camera ID %.0f and IOP [xp yp c] = [%.5f %.5f %.5f] px\n',camID,IOP);

ext = double(S.ext);
imageIDs = ext(:,1);
C = double(S.eop.XYZc);
wpk = double(S.eop.wpk);
assert(size(C,1)==numel(imageIDs) && size(wpk,1)==numel(imageIDs), ...
    'Adjusted EOP row count differs from source EXT image count');
assert(numel(unique(imageIDs))==numel(imageIDs),'Duplicate image IDs');

%% Select original photographs, explicitly by image ID
images = cell(numel(imageIDs),1);
imgWidth = zeros(numel(imageIDs),1);
imgHeight = zeros(numel(imageIDs),1);
imagePaths = strings(numel(imageIDs),1);
for k = 1:numel(imageIDs)
    [ifn,ifd] = uigetfile({'*.jpg;*.jpeg;*.png;*.tif;*.tiff;*.JPG;*.JPEG;*.PNG','Image files'}, ...
       sprintf('Select ORIGINAL photograph for IMAGE %d',imageIDs(k)));
    if isequal(ifn,0),fprintf('Cancelled.\n');return;end
    imagePaths(k) = string(fullfile(ifd,ifn));
    im = imread(char(imagePaths(k)));
    assert(ndims(im)==3 && size(im,3)>=3,'Expected RGB image for image %d',imageIDs(k));
    if ~isa(im,'uint8'), im=im2uint8(im(:,:,1:3)); end
    images{k} = im(:,:,1:3);
    imgHeight(k)=size(im,1); imgWidth(k)=size(im,2);
    fprintf('Image %d: %s (%d x %d)\n',imageIDs(k),ifn,imgWidth(k),imgHeight(k));
end

%% Check camera projection convention against MEASURED IMAGE POINTS
% Use adjusted XYZ for BOTH tie and control points. The control values can
% move slightly because Mode 0 controls are weighted observations.
pointIDs=[double(S.tieRows(:,1));double(S.conRows(:,1))];
pointXYZ=[double(S.tie.XYZ);double(S.con.XYZ)];
assert(numel(unique(pointIDs))==numel(pointIDs), 'Duplicate tie/control IDs');
pho=double(S.pho);
[found,row]=ismember(pho(:,1),pointIDs);
assert(all(found),'Some PHO points have no adjusted 3D XYZ');
obsXYZ=pointXYZ(row,:);
allPixelErr=[];
depthSigns=zeros(numel(imageIDs),1);
fprintf('\nREPROJECTION CHECK (observed pixels against adjusted XYZ):\n');
for k=1:numel(imageIDs)
    q=pho(:,2)==imageIDs(k);
    assert(any(q),'No measured observations for image %d',imageIDs(k));
    [pix,depth]=projectXYZ(obsXYZ(q,:),C(k,:),wpk(k,:),IOP);
    inFrame=isfinite(pix(:,1)) & isfinite(pix(:,2)) & ...
            pix(:,1)>=1 & pix(:,1)<=imgWidth(k) & ...
            pix(:,2)>=1 & pix(:,2)<=imgHeight(k);
    d=pix-pho(q,3:4);
    radial=hypot(d(:,1),d(:,2));
    rmsRad=sqrt(mean(radial.^2,'omitnan'));
    medRad=median(radial,'omitnan');
    allPixelErr=[allPixelErr;radial]; %#ok<AGROW>
    s=sign(median(depth,'omitnan'));
    if s==0 || ~isfinite(s),error('Cannot determine depth sign for image %d',imageIDs(k));end
    depthSigns(k)=s;
    fprintf(' Image %d: median %.3f px, radial RMS %.3f px, inside %d/%d\n', ...
        imageIDs(k),medRad,rmsRad,nnz(inFrame),nnz(q));
end
overallRad=sqrt(mean(allPixelErr.^2,'omitnan'));
fprintf(' Overall radial reprojection RMS: %.3f px\n',overallRad);
if overallRad>5
    error(['Projection check failed: RMS %.2f px exceeds 5 px. ' ...
        'Do not colourize until image/point mapping, rotation or IOPs are corrected.'],overallRad);
end

%% Cropped TLS XYZ input (original scanner datum)
[pfn,pfd]=uigetfile({'*.txt;*.xyz;*.asc;*.csv','Point cloud text (*.txt,*.xyz,*.asc,*.csv)'}, ...
    'Choose cropped CloudCompare point cloud (XYZ or XYZRGB)');
if isequal(pfn,0),return;end
cloudPath=fullfile(pfd,pfn);
answer=menu('Units in cropped TLS point cloud?', ...
    'Metres (most CloudCompare exports)', 'Millimetres');
if answer==0,return;end
if answer==1, mmPerCloudUnit=1000;else,mmPerCloudUnit=1;end

[outfn,outfd]=uiputfile('*.txt','Save coloured XYZRGB point cloud', ...
    fullfile(pfd,'lab5_TLS_colourized_XYZRGB.txt'));
if isequal(outfn,0),return;end
outPath=fullfile(outfd,outfn);
logPath=fullfile(outfd,'lab5_backprojection_report.txt');

% Output colour-selection strategy: among projections within each image,
% choose the smallest normalized distance to the image centre.
% It DOES NOT distinguish front-facing vs occluded surfaces.
inputFID=fopen(cloudPath,'rt');assert(inputFID>0,'Cannot open point cloud');
inputCleanup=onCleanup(@()fclose(inputFID));
outputFID=fopen(outPath,'wt');assert(outputFID>0,'Cannot write output cloud');
outputCleanup=onCleanup(@()fclose(outputFID));
chunkSize=30000;
processed=0;coloured=0;skipped=0;
chosenCount=zeros(numel(imageIDs),1);
lineCount=0;

fprintf('\nProcessing cloud in chunks of %d points.\n',chunkSize);
while true
    lines=cell(chunkSize,1);nRead=0;
    while nRead<chunkSize
        ln=fgetl(inputFID);
        if ~ischar(ln),break;end
        lineCount=lineCount+1;
        ln=strtrim(ln);
        if isempty(ln) || startsWith(ln,'//') || startsWith(ln,'#'),continue;end
        nRead=nRead+1;lines{nRead}=ln;
    end
    if nRead==0,break;end
    X=zeros(nRead,3); valid=false(nRead,1);
    for t=1:nRead
        % Accept whitespace, comma or semicolon separators; ignore extra RGB fields.
        nums=sscanf(regexprep(lines{t},'[,;]',' '),'%f');
        if numel(nums)>=3 && all(isfinite(nums(1:3)))
            X(t,:)=nums(1:3)';valid(t)=true;
        elseif processed==0 && t<=5
            fprintf('Skipping non-numeric/header line: %s\n',lines{t});
        end
    end
    X=X(valid,:); if isempty(X),continue;end
    N=size(X,1);
    XYZmm=X*mmPerCloudUnit;
    RGB=zeros(N,3,'uint8');
    bestScore=inf(N,1);
    bestCamera=zeros(N,1);
    for k=1:numel(imageIDs)
        [pix,depth]=projectXYZ(XYZmm,C(k,:),wpk(k,:),IOP);
        col=round(pix(:,1));row=round(pix(:,2));
        validProj=isfinite(pix(:,1)) & isfinite(pix(:,2)) & ...
          (depth*depthSigns(k)>0) & ...
          col>=1 & col<=imgWidth(k) & row>=1 & row<=imgHeight(k);
        if ~any(validProj),continue;end
        % Normalized center score allows photographs of slightly varying sizes.
        score=((pix(:,1)-(imgWidth(k)+1)/2)/(imgWidth(k)/2)).^2 + ...
              ((pix(:,2)-(imgHeight(k)+1)/2)/(imgHeight(k)/2)).^2;
        take=find(validProj & score<bestScore);
        if isempty(take),continue;end
        linearPixel=sub2ind([imgHeight(k),imgWidth(k)],row(take),col(take));
        im=images{k};
        for channel=1:3
            a=im(:,:,channel);
            RGB(take,channel)=a(linearPixel);
        end
        bestScore(take)=score(take);
        bestCamera(take)=k;
    end
    keep=bestCamera>0;
    nn=nnz(keep);coloured=coloured+nn;processed=processed+N;skipped=skipped+(N-nn);
    for k=1:numel(imageIDs),chosenCount(k)=chosenCount(k)+nnz(bestCamera==k);end
    if nn>0
        % Write exactly six columns, retaining original point cloud units.
        out=[X(keep,:),double(RGB(keep,:))];
        fprintf(outputFID,'%.8f %.8f %.8f %d %d %d\n',out');
    end
    if mod(processed,300000)<chunkSize
        fprintf(' Processed %d, coloured %d, outside all images %d\n',processed,coloured,skipped);
    end
end

%% Summary log
lf=fopen(logPath,'wt');
if lf>0
    fprintf(lf,'ENGO 531 Lab 5 - Backprojection / nearest neighbour RGB\n');
    fprintf(lf,'Input MAT: %s\nInput cloud: %s\nOutput cloud: %s\n', ...
        fullfile(fd,fn),cloudPath,outPath);
    fprintf(lf,'IOP xp yp c (px): %.7f %.7f %.7f\n',IOP);
    fprintf(lf,'Projection rotation: (Rx*Ry*Rz) transpose; x=xp-c*U/W; y=yp+c*V/W\n');
    fprintf(lf,'Pixel interpolation: nearest neighbor\n');
    fprintf(lf,'Selection: projection closest to normalized image centre\n');
    fprintf(lf,'No occlusion test; points inside image bounds can still be hidden\n');
    fprintf(lf,'Reprojection check radial RMS %.4f px\n',overallRad);
    fprintf(lf,'Original cloud units: %s\n',stringUnit(mmPerCloudUnit));
    fprintf(lf,'Total numeric input points: %d\nColoured points: %d\nExcluded points: %d\n', ...
        processed,coloured,skipped);
    for k=1:numel(imageIDs)
        fprintf(lf,'Image %d assigned points %d: %s\n',imageIDs(k),chosenCount(k),imagePaths(k));
    end
    fclose(lf);
end
fprintf('\nFinished. Coloured %d/%d valid points (%.1f%%).\n', ...
    coloured,processed,100*coloured/max(1,processed));
fprintf('CloudCompare six-column XYZRGB: %s\n',outPath);
fprintf('Processing report: %s\n',logPath);
fprintf('IMPORTANT: visually inspect seams and wrong colours: image-center selection does not test occlusion.\n');
end

function [pixel,depth]=projectXYZ(XYZ,cameraXYZ,wpk,IOP)
    a=deg2rad(wpk(1));b=deg2rad(wpk(2));g=deg2rad(wpk(3));
    Rx=[1 0 0;0 cos(a) -sin(a);0 sin(a) cos(a)];
    Ry=[cos(b) 0 sin(b);0 1 0;-sin(b) 0 cos(b)];
    Rz=[cos(g) -sin(g) 0;sin(g) cos(g) 0;0 0 1];
    M=(Rx*Ry*Rz)';
    q=(M*(XYZ-cameraXYZ)')';
    depth=q(:,3);
    d=depth; d(abs(d)<1e-10)=NaN;
    pixel=[IOP(1)-IOP(3)*q(:,1)./d, IOP(2)+IOP(3)*q(:,2)./d];
end

function s=stringUnit(f)
    if f==1000,s='metres';else,s='millimetres';end
end
