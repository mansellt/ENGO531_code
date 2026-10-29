%% ENGO 531 Lab 4 - FIELD NOTE EOPs, reversed phi, images 1-5 (minimum 2)
% Uses documented ENGO 531 Lab 1 file formats and operation modes.
% Parameter observations constrain ONLY the FIRST IMAGE EOPs (Lab 1 p3).
% Camera IOPs remain fixed. No data snooping.
% Camera positions, object coordinates and distance are in mm.
clear; clc; close all;

%% USER SETTINGS
imageWidth = 5712; imageHeight = 4284;
IOP = [2834.0 2156.7 3995.9];   % calibrated Lab 3 xp, yp, c [pixels]
sd_xy = 3.2;                     % assumed image coordinate SD [pixels]
sd_dist = 5;                     % tentative distance SD [mm]; replace with tape SD
sd_firstEOP_XYZ = 1;             % mm, tight datum on first image EOP position
sd_firstEOP_wpk = 0.01;          % deg, tight datum on first image EOP angles
scaleFrom = "F1"; scaleTo = "F6";
measuredDistance_mm = 2478;       % 2.478 metres
minImagesPerPoint = ;
allowF1TwoImageTrial = true;      % F1 meets the 2-image minimum
runBundle = true;
trialImageIDs = [1 2 3 4 5 6 7 8 9]; % connected component; images 7-9 have no common points

% Approximate EOPs from your field notes:
% Columns: image, X, Y, Z [mm], phi, kappa, omega [degrees]
notes = [
 1 -2213  1711 -645  -90  0 90;
 2 -2279   561 -495  -75  0 90;
 3 -1619 -1540 -495  -45  0 90;
 4  -440 -2027 -495  -25  0 90;
 5  1239 -2545 -495    0  0 90;
 6  3486 -1652 -495   30  0 90;
 7  4368   -97 -495   60  0 90;
 8  4901   411 -495   90  0 90;
 9  4455  2686 -495  120  0 90];

%% SELECT FILES
[fn,fd] = uigetfile({'*.PHO;*.pho','Measured image points'}, ...
    'Select calibration(2).PHO');
if isequal(fn,0), return; end
phoFile = fullfile(fd,fn);
outputDir = uigetdir(fd,'Select folder for Lab 4 trial output');
if isequal(outputDir,0), return; end
prefix = fullfile(outputDir,'engo531_lab4_fieldEOP_negativePhi_min2');

% PHO point IDs may be alphanumeric; Lab 1 explicitly supports them.
% We map to numeric IDs for stable file generation and report the map.
fid = fopen(phoFile,'r'); assert(fid>=0,'Could not open PHO');
C = textscan(fid,'%s %f %f %f'); fclose(fid);
pointName = string(C{1}); imgID = C{2}; xy = [C{3},C{4}];
assert(~isempty(imgID) && all(isfinite([imgID;xy(:)])),'Invalid PHO');
assert(all(xy(:,1)>=0 & xy(:,1)<=imageWidth & ...
           xy(:,2)>=0 & xy(:,2)<=imageHeight),'Points outside image bounds');
assert(numel(unique(pointName+"_"+string(imgID)))==numel(imgID), ...
       'Duplicate point-image observations');

% Image 6 is in field notes but not the supplied PHO; use actual image set.
activeImages = intersect(notes(:,1),unique(imgID),'stable');
activeImages = activeImages(ismember(activeImages,trialImageIDs));
missingImages = setdiff(notes(:,1),activeImages);
fprintf('Images used: %s\n',mat2str(activeImages'));
if ~isempty(missingImages)
    warning('Images without PHO measurements excluded: %s',mat2str(missingImages'));
end

% Filter PHO to selected connected image block BEFORE counting overlap.
blockSel = ismember(imgID,activeImages);
pointName=pointName(blockSel); imgID=imgID(blockSel); xy=xy(blockSel,:);

% Retain all points appearing in TWO or more selected images.
[names,~,group] = unique(pointName,'stable');
counts = accumarray(group,1);
disp(table(names,counts,'VariableNames',{'Target','ImageCount'}));
keepNames = names(counts>=minImagesPerPoint);
endpoints = [scaleFrom;scaleTo];
for t=1:2
    id = find(names==endpoints(t),1);
    assert(~isempty(id),'Scale endpoint %s not measured in PHO',endpoints(t));
    if ~ismember(endpoints(t),keepNames)
        if endpoints(t)=="F1" && counts(id)==2 && allowF1TwoImageTrial
            warning('F1 has two observations and is retained under the requested minimum');
            keepNames(end+1)=endpoints(t); %#ok<SAGROW>
        else
            error('Scale point %s has %d images; minimum %d required',endpoints(t),counts(id),minImagesPerPoint);
        end
    end
end
sel=ismember(pointName,keepNames) & ismember(imgID,activeImages);
pointName=pointName(sel); imgID=imgID(sel); xy=xy(sel,:);
[names,~,ptIndex]=unique(pointName,'stable');
[tf,imgIndex]=ismember(imgID,activeImages); assert(all(tf));
assert(all(ismember(endpoints,names)),'Scale endpoints absent after filtering');
npts=numel(names); nimg=numel(activeImages); nobs=numel(imgID);
fprintf('Retained: %d image points, %d targets, %d images\n',nobs,npts,nimg);

% Persist reversible numeric IDs for the course file interface.
% F1 -> 1, F6 -> 6 etc are NOT assumed; numeric IDs are assigned explicitly.
numericID = (1001:1000+npts)';
idMap=table(names,numericID,'VariableNames',{'OriginalPointID','BundlePointID'});
writetable(idMap,[prefix '_point_id_map.csv']);
[tf,scaleIdx]=ismember(endpoints,names); assert(all(tf));

% Convert phi/kappa/omega field-note order to omega/phi/kappa course order.
[tf,r]=ismember(activeImages,notes(:,1)); assert(all(tf));
N=notes(r,:);
ext=[N(:,1),repmat(99,nimg,1),N(:,2:4),N(:,7),N(:,5),N(:,6)];

% INITIAL TIE POINT XYZ FROM PREVIOUS JOINT IMAGE/OBJECT PRE-FIT
% Previous numerical least-squares initialization minimized all
% 40 measured image-point reprojection discrepancies simultaneously,
% keeping Image 1 camera pose fixed and constraining F1--F6 to 2478 mm.
% These starting values are approximate; NOT final adjusted coordinates.
% No image observations were rejected or snooped.
% The P-code still performs the course adjustment on the original PHO.
% Reprojection diagnostic uses a CANDIDATE rotation model.
% Its equivalence to the P-code for arbitrary angles is not established.
initialXYZNames = [
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
    "L9";
];
initialXYZ = [
    -50.9988245 48.8389883 -38.7036443;
    -86.4568735 51.3483000 -2133.2035782;
    702.9923490 46.4639473 -51.3504155;
    1674.0877894 75.8147053 -70.2513596;
    2425.8978800 116.1343280 -81.8228151;
    2418.7110316 149.9278108 -2181.5225131;
    -153.2735849 591.1482736 -1078.9046793;
    -173.5873315 735.2712955 -1781.2878620;
    -192.6145164 1821.9242191 -37.6329990;
    -114.3387807 16.0720985 -57.5115687;
    -147.9520702 33.0009521 -2117.0338348;
    -225.3004422 1839.6915304 -2101.7075262;
    -186.6556773 1412.0734304 -1071.6515515;
];
% Known 13 joint-pre-fit coordinates are retained. Additional two-image
% targets get geometric approximations from the labeled image planes.
% Geometry is for INITIAL VALUES only and is NOT treated as fixed control.
[ok, idxXYZ] = ismember(names,initialXYZNames);
XYZ0 = zeros(npts,3);
XYZ0(ok,:) = initialXYZ(idxXYZ(ok),:);
if any(~ok)
    fid=fopen(phoFile,'r');
    raw=textscan(fid,'%s %f %f %f'); fclose(fid);
    allName=string(raw{1}); allImg=raw{2}; allXY=[raw{3},raw{4}];
    Hleft=planeMap(["L1" "L4" "L5" "L8"],1,allName,allImg,allXY, ...
         [0 0;2704 0;2704 -2213;0 -2213]);
    Hfront4=planeMap(["F1" "F6" "F7" "F12"],4,allName,allImg,allXY, ...
         [0 0;2478 0;2478 -2213;0 -2213]);
    Hfront5=planeMap(["F1" "F6" "F7" "F12"],5,allName,allImg,allXY, ...
         [0 0;2478 0;2478 -2213;0 -2213]);
    for j=find(~ok)'
        nm=names(j);
        if startsWith(nm,"L")
            q=find(allName==nm & allImg==1,1);
            if isempty(q) && nm=="L12"
                XYZ0(j,:)=[-65 1352 -2213];
                warning('L12 absent from image 1: using approximate XYZ');
            else
                assert(~isempty(q),'No Image 1 observation for %s',nm);
                uv=mapPixel(Hleft,allXY(q,:));
                XYZ0(j,:)=[-65,2704-uv(1),uv(2)];
            end
        elseif startsWith(nm,"F")
            q=find(allName==nm & allImg==5,1);
            H=Hfront5;
            if isempty(q)
                q=find(allName==nm & allImg==4,1);
                H=Hfront4;
            end
            assert(~isempty(q),'No image 4/5 observation for %s',nm);
            uv=mapPixel(H,allXY(q,:));
            XYZ0(j,:)=[uv(1),0,uv(2)];
        else
            error('Target %s missing usable initial XYZ',nm);
        end
    end
end
fprintf('Initial XYZ values: %d joint-prefit and %d geometric additions.\n', ...
    sum(ok),sum(~ok));

% EOPs are taken DIRECTLY from the user's field-note table above.
% The .ext columns are [imageID cameraID Xc Yc Zc omega phi kappa].
fprintf('\nField-note camera XYZ with reversed phi angles:\n');
disp(array2table(ext(:,[1,3:8]),'VariableNames', ...
 {'Image','Xc_mm','Yc_mm','Zc_mm','omega_deg','phi_deg','kappa_deg'}));
fprintf('Initial F1-F6 scale distance = %.3f mm\n', ...
 norm(XYZ0(scaleIdx(1),:)-XYZ0(scaleIdx(2),:)));

%% PRE-ADJUSTMENT REPROJECTION DIAGNOSTICS
% Diagnostic projection only: DOES NOT establish P-code convention.
% The prefit XYZ were fitted with different camera rotations; large initial
% reprojection discrepancies are expected under this intentional test.
errors = zeros(nobs,1);
for h=1:nobs
    e=ext(imgIndex(h),:);
    M=rotationWPK_approx(deg2rad(e(6:8)));
    pc=M*(XYZ0(ptIndex(h),:)'-e(3:5)');
    if abs(pc(3))<1e-9
        errors(h)=Inf;
    else
        predicted=[IOP(1)-IOP(3)*pc(1)/pc(3); IOP(2)+IOP(3)*pc(2)/pc(3)];
        errors(h)=norm(predicted-xy(h,:)');
    end
end
fprintf('Initialization: median reprojection discrepancy %.1f px; RMS %.1f px\n', ...
    median(errors),sqrt(mean(errors.^2)));
for k=1:nimg
    ss=imgIndex==k;
    fprintf('  Image %d median reprojection discrepancy: %.1f px\n', ...
        activeImages(k),median(errors(ss)));
end
if median(errors)>100
    warning('Initial discrepancy exceeds 100 px; P-code convergence is NOT guaranteed.');
end

% No control coordinates: datum via EOP parameter observations, scale via DIS.
% Both point types are free XYZ, stored as .tie.
phoOut=[numericID(ptIndex),imgID,xy];
tieOut=[numericID,XYZ0];
writeNumeric([prefix '.pho'],phoOut,'%.0f %.0f %.2f %.2f\n');
writematrix(ext,[prefix '.ext'],'Delimiter','tab','FileType','text');
writematrix(tieOut,[prefix '.tie'],'Delimiter','tab','FileType','text');

% Same INT format used by the successful Lab 3 course script:
% camera ID, xp, yp, c, -1 (no distortion terms).
int=[99 IOP -1];
writematrix(int,[prefix '.int'],'Delimiter','tab','FileType','text');

% Lab 1 p.3: CFG nine lines: param_obs_mode, iop_mode,
% tolerances XYZ/wpk/IOP, SD image/distance/EOP XYZ/EOP angles.
% param_obs_mode=1 -> parameter observations for SIX EOPs of FIRST image ONLY.
% iop_mode=0 -> IOPs are fixed constants (not estimated).
cfg=[1;0;0.001;0.0001;0.01;sd_xy;sd_dist; ...
     sd_firstEOP_XYZ;sd_firstEOP_wpk];
writematrix(cfg,[prefix '.cfg'],'FileType','text');

% Lab 1 p.2: DIS has three columns: from_point_id to_point_id distance.
% Distance uses SAME units as EXT and TIE coordinates (here, mm).
dis=[numericID(scaleIdx(1)),numericID(scaleIdx(2)),measuredDistance_mm];
writeNumeric([prefix '.dis'],dis,'%.0f %.0f %.4f\n');
if isfile([prefix '.con']),delete([prefix '.con']);end

% Diagnostics that help assess whether the network is solvable.
unknowns=6*nimg+3*npts;
observations=2*nobs+1+6; % 6 first-image EOP observations, 1 scale distance
dof=observations-unknowns;
fprintf('Unknowns = %d; scalar observations = %d; nominal DOF = %d\n', ...
    unknowns,observations,dof);
assert(dof>0,'Nonpositive network redundancy after filtering');
fprintf('First active image %d gets six weighted EOP observations.\n',activeImages(1));
fprintf('Input prefix: %s\n',prefix);
fprintf('Using fixed IOP [xp yp c]=[%.3f %.3f %.3f] px\n',IOP);
fprintf('NO data snooping or automatic residual deletion.\n');
fprintf('Tie-point XYZ from previous joint prefit; field-note camera XYZ and reversed phi angles.\n');
fprintf('EOPs match requested original XYZ, omega=90, kappa=0, phi=-90,-75,-45,-25,0.\n');

if runBundle
    assert(exist('full_matlab_bundle','file')==6 || ...
           exist('full_matlab_bundle','file')==2, ...
      'Rename full_matlab_bundle(3).p to full_matlab_bundle.p and add its folder to MATLAB path');
    try
        [eop,iop,conAdjusted,tie,Cx,vr_img,vr_eop,vr_dis,vr_con] = ...
            full_matlab_bundle(prefix);
    catch ME
        fprintf(2,'Course bundle failed: %s\n',ME.message);
        fprintf(2,'Inputs follow Lab 1 specifications. Check initial XYZ rays and angle conventions.\n');
        rethrow(ME);
    end
    fprintf('\nBundle returned (check MATLAB convergence message). Image residual RMS: %.4f pixels\n', ...
        sqrt(mean(vr_img(:,1).^2)));
    disp('Adjusted EOPs [mm, degrees]:');disp(eop.XYZc);disp(eop.wpk);
    if isfield(tie,'XYZ'),disp('Adjusted target XYZ [mm]:');disp(tie.XYZ);end
    if ~isempty(vr_dis),fprintf('Distance residual (mm): %.4f\n',vr_dis(1,1));end
    if ~isempty(vr_eop),fprintf('First-image EOP residuals:\n');disp(vr_eop(:,1));end
    save([prefix '_results.mat'],'eop','iop','conAdjusted','tie','Cx', ...
        'vr_img','vr_eop','vr_dis','vr_con','idMap','ext','int','cfg', ...
        'XYZ0','phoOut','tieOut');
end

function R=rotationWPK_approx(wpk)
w=wpk(1);p=wpk(2);k=wpk(3);
Rx=[1 0 0;0 cos(w) -sin(w);0 sin(w) cos(w)];
Ry=[cos(p) 0 sin(p);0 1 0;-sin(p) 0 cos(p)];
Rz=[cos(k) -sin(k) 0;sin(k) cos(k) 0;0 0 1];
R=Rz'*Ry*Rx'; % candidate course model; confirm M with course software
end
function writeNumeric(path,A,fmt)
fid=fopen(path,'w');assert(fid>=0,'Could not open %s',path);
c=onCleanup(@()fclose(fid));
for i=1:size(A,1),fprintf(fid,fmt,A(i,:));end
end


function H=planeMap(corners,imageID,allNames,allImages,allXY,uv)
A=zeros(8,8);b=zeros(8,1);
for ii=1:4
    q=find(allNames==corners(ii) & allImages==imageID,1);
    assert(~isempty(q),'Anchor %s missing in image %d',corners(ii),imageID);
    x=allXY(q,1);y=allXY(q,2);u=uv(ii,1);v=uv(ii,2);
    A(2*ii-1,:)=[x y 1 0 0 0 -u*x -u*y];
    A(2*ii,:)=[0 0 0 x y 1 -v*x -v*y];
    b(2*ii-1:2*ii)=[u;v];
end
assert(rcond(A)>1e-14,'Degenerate planar map');
h=A\b;
H=[h(1) h(2) h(3);h(4) h(5) h(6);h(7) h(8) 1];
end
function uv=mapPixel(H,xy)
p=H*[xy(:);1];
assert(abs(p(3))>1e-10,'Bad planar mapping denominator');
uv=(p(1:2)/p(3))';
end
