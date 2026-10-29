%% ENGO 531 Lab 4 - course full_matlab_bundle.p trial (Lab 1 verified formats)
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
minImagesPerPoint = 3;
allowF1TwoImageTrial = true;      % exception ONLY to get a preliminary run
runBundle = true;
trialImageIDs = [1 2 3 4 5]; % connected component; images 7-9 have no common points

% Approximate EOPs from your field notes:
% Columns: image, X, Y, Z [mm], phi, kappa, omega [degrees]
notes = [
 1 -2213  1711 -645   90    0 90;
 2 -2279   561 -495   75    0 90;
 3 -1619 -1540 -495   45    0 90;
 4  -440 -2027 -495   25    0 90;
 5  1239 -2545 -495    0    0 90;
 6  3486 -1652 -495  -30    0 90;
 7  4368   -97 -495  -60    0 90;
 8  4901   411 -495  -90    0 90;
 9  4455  2686 -495 -120    0 90];

%% SELECT FILES
[fn,fd] = uigetfile({'*.PHO;*.pho','Measured image points'}, ...
    'Select calibration(2).PHO');
if isequal(fn,0), return; end
phoFile = fullfile(fd,fn);
outputDir = uigetdir(fd,'Select folder for Lab 4 trial output');
if isequal(outputDir,0), return; end
prefix = fullfile(outputDir,'engo531_lab4_course_trial_verified');

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

% Minimum three image observations per target, with optional trial exception.
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
            warning('TRIAL EXCEPTION: F1 has only two observations');
            keepNames(end+1)=endpoints(t); %#ok<SAGROW>
        else
            error('Scale point %s has %d images; at least 3 required',endpoints(t),counts(id));
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

% Approximate 3D points: linear intersection of measured image rays.
% These are STARTING VALUES, not fixed XYZ observations.
% The rotation convention MUST match the one used by full_matlab_bundle.p.
XYZ0=zeros(npts,3);
for j=1:npts
    ii=find(ptIndex==j);
    A=zeros(3); b=zeros(3,1);
    for h=ii'
        e=ext(imgIndex(h),:);
        Cc=e(3:5)';
        R=rotationWPK_approx(deg2rad(e(6:8)));
        rayCam=[xy(h,1)-IOP(1); xy(h,2)-IOP(2); -IOP(3)];
        rayWorld=R'*rayCam; rayWorld=rayWorld/norm(rayWorld);
        M=eye(3)-rayWorld*rayWorld';
        A=A+M; b=b+M*Cc;
    end
    assert(rcond(A)>1e-12,'Poor ray intersection for %s',names(j));
    XYZ0(j,:)=(A\b)';
end
fprintf('Approx. F1-F6 distance = %.3f mm; measured = %.3f mm\n', ...
    norm(XYZ0(scaleIdx(1),:)-XYZ0(scaleIdx(2),:)),measuredDistance_mm);

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
fprintf('Ray initialization uses M=Rx(-omega)*Ry(-phi)*Rz(-kappa), a candidate matching Lab 3 reference at 90,0,0.\n');
fprintf('full course rotation order remains to be verified if convergence fails.\n');

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
R=Rx'*Ry'*Rz'; % MATCHES course sample M at w=90,p=k=0; verify full convention
end
function writeNumeric(path,A,fmt)
fid=fopen(path,'w');assert(fid>=0,'Could not open %s',path);
c=onCleanup(@()fclose(fid));
for i=1:size(A,1),fprintf(fid,fmt,A(i,:));end
end
