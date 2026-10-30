function ENGO531_Lab4_combine_1to9()
% ENGO 531 LAB 4 — ONE JOINT 9-IMAGE BUNDLE, initialized from TWO
% CONVERGED independent course adjustments (images 1-6 and images 7-9).
%
% Requires:
%   (1) *_results.mat from converged images 1-6
%   (2) *_results.mat from converged right-wall images 7-9
%   (3) ORIGINAL all-nine-image .tie with matching .pho/.ext/.int/.cfg/.dis,
%       and its *_target_ID_map.csv, made by ENGO531_Lab4_select_images
%       with selectedImages=1:9, minImagesPerPoint=2 (or 3).
%   (4) Course full_matlab_bundle.p on MATLAB path.
%
% Geometry: resect Image 6 against the converged right-wall model,
% then register the RIGHT frame to the LEFT frame via Image 6's pose.
% Both camera positions AND camera rotations are transformed.
% This alters STARTING VALUES ONLY. All original PHO observations are kept.
%
% Caution: Planar resection can have alternate minima. If the preflight
% checks fail, the joint adjustment will NOT run.
clc;
maxResectionRMS = 20; % pixels; loose preflight threshold for matching
maxJointRMS = 150;    % pixels; allows nonlinear refinement, but blocks gross mistakes
runAnyway = false;    % set true ONLY to experiment despite failed preflight

[fn,fd]=uigetfile('*_results.mat','Select CONVERGED Images 1-6 results');
if isequal(fn,0),return;end
A=load(fullfile(fd,fn));
[fn,fd]=uigetfile('*_results.mat','Select CONVERGED Images 7-9 results');
if isequal(fn,0),return;end
B=load(fullfile(fd,fn));
[fn,fd]=uigetfile('*.tie','Select ORIGINAL all-nine-image .tie file');
if isequal(fn,0),return;end
srcPrefix=fullfile(fd,fn(1:end-4));
[outDirOK,outDir]=getOutputDirectory(fd);
if ~outDirOK,return;end
outPrefix=fullfile(outDir,'engo531_all9_two_network_warmstart');

assert(isfield(A,'eop') && isfield(A,'tie') && isfield(A,'names') && isfield(A,'numericID'), ...
    'Left network MAT does not contain expected fields.');
assert(isfield(B,'eop') && isfield(B,'tie') && isfield(B,'names'), ...
    'Right network MAT does not contain expected fields.');
assert(size(A.eop.XYZc,1)==6 && size(B.eop.XYZc,1)==3, ...
    'Select the SIX-image and THREE-image successful results, in that order.');
assert(isfield(A.tie,'XYZ') && isfield(B.tie,'XYZ'),'Missing adjusted tie.XYZ');

exts={'.pho','.ext','.tie','.int','.dis','.cfg'};
for k=1:numel(exts)
    assert(isfile([srcPrefix exts{k}]),'Missing original file: %s',[srcPrefix exts{k}]);
end
mapFile=[srcPrefix '_target_ID_map.csv'];
assert(isfile(mapFile),'Missing ID map for original all-nine input');
map=readtable(mapFile,'TextType','string');
assert(all(ismember({'OriginalPointID','BundlePointID'},map.Properties.VariableNames)), ...
    'Target ID map does not have expected column names.');
T=readmatrix([srcPrefix '.tie'],'FileType','text');
E=readmatrix([srcPrefix '.ext'],'FileType','text');
P=readmatrix([srcPrefix '.pho'],'FileType','text');
assert(size(T,2)==4 && size(E,2)==8 && size(P,2)==4,'Unexpected course file columns');
assert(isequal(E(:,1),(1:9)'),'The all-nine .ext must contain images 1 through 9, in order.');
assert(numel(unique(E(:,1)))==9,'Repeated camera records.');
[tf,tmap]=ismember(T(:,1),double(map.BundlePointID));assert(all(tf));
labels=string(map.OriginalPointID(tmap));

% The saved tie rows correspond to input .tie IDs; verify where available.
AL=string(A.names(:)); BL=string(B.names(:));
assert(numel(AL)==size(A.tie.XYZ,1) && numel(BL)==size(B.tie.XYZ,1), ...
    'Returned adjusted ties and target labels do not align by length.');
assertCourseIDs(A,'numericID');
assertCourseIDs(B,'id');

% Only right-wall targets in the right model participate in resection.
[tfB,bi]=ismember(labels,BL);
Rmask=startsWith(labels,'R') & tfB;
assert(nnz(Rmask)>=8,'Too few common R targets between all-nine file and right-wall MAT.');
% Obtain Image-6 original image observations, selecting R points only.
rows=(P(:,2)==6) & ismember(P(:,1),T(Rmask,1));
img6=P(rows,:); assert(size(img6,1)>=8,'Need at least 8 Image-6 R observations for resection.');
[tf,rt]=ismember(img6(:,1),T(:,1)); assert(all(tf));
[tf,rr]=ismember(labels(rt),BL);assert(all(tf));
Xright6=B.tie.XYZ(rr,:);
obs6=img6(:,3:4);
% Camera convention verified previously against adjusted 1-6 solution:
% M = (Rz(kappa)*Ry(phi)*Rx(omega))'', depth = camera Z;
% x = xp - c*Qx/Qz; y = yp + c*Qy/Qz.
IOP=[2834.0 2156.7 3995.9];
if isfield(A,'IOP'),IOP=double(A.IOP(:)');end
if isfield(B,'IOP')
 assert(max(abs(double(B.IOP(:)')-IOP))<1e-4, ...
    'The two solutions used different IOPs; use the same calibrated IOPs.');
end

% Check our camera convention reproduces the successful six-image fit.
% Build a test from the six-image MAT's own saved photoTable if available.
if isfield(A,'photoTable')
 q=A.photoTable; idA=A.numericID(:);
 [has,ix]=ismember(q(:,1),idA);
 assert(all(has));
 pred=zeros(size(q,1),2);
 for i=1:6
  rowsA=q(:,2)==i;
  pose=[A.eop.XYZc(i,:),A.eop.wpk(i,:)];
  pred(rowsA,:)=project(A.tie.XYZ(ix(rowsA),:),pose,IOP);
 end
 refRMS=sqrt(mean(sum((pred-q(:,3:4)).^2,2)));
 fprintf('CHECK: converged 1-6 reference reprojection RMS (radial) = %.3f px\n',refRMS);
 assert(isfinite(refRMS)&&refRMS<15,'Projection convention or tie ordering validation FAILED.');
else
 fprintf('Reference photoTable absent; skipping independent projection-convention test.\n');
end

% Resection image6 in right-wall datum. First guess from historical EOP.
initials=[3941.8 -933 170 67.23 44.50 26.37; ...
          3486 -1652 -495 90 30 0; ...
          4368 -97 -495 90 60 0; ...
          3600 -1200 0 90 45 0];
[bestPose,bestRMS,trialLog]=solvePose(Xright6,obs6,IOP,initials);
fprintf('\nIMAGE 6 RESECTION IN RIGHT-WALL REFERENCE FRAME\n');
disp(array2table(trialLog,'VariableNames',{'Seed','RMS_Pixels','Iterations'}));
fprintf('Image 6 pose: XYZ [%.2f %.2f %.2f] mm; wpk [%.3f %.3f %.3f] deg\n',bestPose);
fprintf('Resection radial RMS %.3f px\n',bestRMS);
dxy=project(Xright6,bestPose,IOP)-obs6;
errors=hypot(dxy(:,1),dxy(:,2));
writetable(sortrows(table(labels(rt),errors,'VariableNames',{'Target','ResectionRadialError_px'}), ...
 'ResectionRadialError_px','descend'),[outPrefix '_image6_resection_errors.csv']);
if bestRMS>maxResectionRMS && ~runAnyway
 fprintf(2,'STOP: Image 6 resection RMS %.1f > %.1f px. Alignment is not reliable.\n', ...
    bestRMS,maxResectionRMS);
 fprintf('Inspect %s_image6_resection_errors.csv before proceeding.\n',outPrefix);
 return;
end

% Relationship Xleft = R * Xright + t, using Image 6's matched poses.
pose6Left=[A.eop.XYZc(6,:),A.eop.wpk(6,:)];
ML=rotation(pose6Left(4:6)); MR=rotation(bestPose(4:6));
R=ML'*MR; t=pose6Left(1:3)'-R*bestPose(1:3)';
fprintf('\nFRAME ALIGNMENT (RIGHT -> LEFT)\n');
fprintf('det(R) %.8f, translation [%.2f %.2f %.2f] mm\n',det(R),t);

Tnew=T; Enew=E;
% L/F use converged original solution. All R use transformed right solution.
reusedLeft=0; reusedRight=0;
for j=1:size(Tnew,1)
 if startsWith(labels(j),'R')
  idx=find(BL==labels(j),1);
  if ~isempty(idx)
   Tnew(j,2:4)=(R*B.tie.XYZ(idx,:)'+t)';
   reusedRight=reusedRight+1;
  end
 else
  idx=find(AL==labels(j),1);
  if ~isempty(idx)
   Tnew(j,2:4)=A.tie.XYZ(idx,:);
   reusedLeft=reusedLeft+1;
  end
 end
end
assert(reusedRight>=12 && reusedLeft>=15,'Too few adjusted tie points reused.');
% Images 1-6 in the left datum as solved.
Enew(1:6,3:5)=A.eop.XYZc;
Enew(1:6,6:8)=A.eop.wpk;
% Images 7-9: both camera centers and orientations MUST be transformed.
for k=7:9
 kR=k-6;
 Cr=B.eop.XYZc(kR,:)';
 Enew(k,3:5)=(R*Cr+t)';
 Mr=rotation(B.eop.wpk(kR,:));
 Ml=Mr*R';
 Enew(k,6:8)=anglesFromRotation(Ml);
end
% Verify transformed right model reproduces original right photo measurements.
[has,pointrow]=ismember(P(:,1),Tnew(:,1));assert(all(has));
report=zeros(9,3);
for k=1:9
 q=P(:,2)==k;
 pred=project(Tnew(pointrow(q),2:4),Enew(k,3:8),IOP);
 residual=pred-P(q,3:4);
 e=hypot(residual(:,1),residual(:,2));
 report(k,:)=[k,sum(q),sqrt(mean(e.^2))];
end
fprintf('\nBEFORE-BUNDLE REPROJECTION VALIDATION\n');
disp(array2table(report,'VariableNames',{'Image','PhotoPoints','RadialRMS_px'}));
writetable(array2table(report,'VariableNames',{'Image','PhotoPoints','RadialRMS_px'}), ...
 [outPrefix '_preflight_per_image.csv']);
if max(report(:,3))>maxJointRMS && ~runAnyway
 fprintf(2,'STOP: Some pre-bundle image RMS exceeds %.1f px; check alignment.\n',maxJointRMS);
 return;
end

% Input file copies ensure original measurements remain unchanged.
for k=1:numel(exts)
 copyfile([srcPrefix exts{k}],[outPrefix exts{k}]);
end
writeRows([outPrefix '.tie'],Tnew,'%d %.7f %.7f %.7f\n');
writeRows([outPrefix '.ext'],Enew,'%d %d %.7f %.7f %.7f %.8f %.8f %.8f\n');
% Keep the ORIGINAL F1--F6 distance and first-camera datum from source.
% The right-wall R4--R1 distance is already reflected in the initialized
% right reconstruction and can be used as an independent validation.
copyfile(mapFile,[outPrefix '_target_ID_map.csv']);
writetable(table(labels,Tnew(:,2),Tnew(:,3),Tnew(:,4), ...
 'VariableNames',{'Target','StartX_mm','StartY_mm','StartZ_mm'}), ...
 [outPrefix '_starting_XYZ.csv']);
fprintf('\nPrepared joint adjustment: %s.*\n',outPrefix);
fprintf('Reused %d adjusted L/F points and %d adjusted R points.\n',reusedLeft,reusedRight);
fprintf('Original PHO was COPIED WITHOUT CHANGES.\n');

reply=input('Run the FULL 9-image course bundle now? [y/N]: ','s');
if ~strcmpi(strtrim(reply),'y'),return;end
assert(exist('full_matlab_bundle','file')==2 || exist('full_matlab_bundle','file')==6, ...
 'full_matlab_bundle.p not on MATLAB path');
[eop,iop,con,tie,Cx,vr_img,vr_eop,vr_dis,vr_con]=full_matlab_bundle(outPrefix);
rmsXY=sqrt(mean(vr_img(:,1).^2));
fprintf('Returned one-component RMS = %.3f px (valid ONLY if course reports convergence)\n',rmsXY);
fprintf('Returned camera XYZ [mm]:\n');disp(eop.XYZc);
fprintf('Returned camera omega/phi/kappa [deg]:\n');disp(eop.wpk);
save([outPrefix '_results.mat'],'eop','iop','con','tie','Cx','vr_img','vr_eop', ...
 'vr_dis','vr_con','R','t','bestPose','bestRMS','Tnew','Enew','labels','report','rmsXY');
fprintf('Verify the course explicitly printed: Adjustment converged.\n');
end

function [ok,out]=getOutputDirectory(defaultDir)
out=uigetdir(defaultDir,'Select folder to SAVE combined 9-image trial');
ok=~isequal(out,0);
end

function assertCourseIDs(S,name)
if ~isfield(S,name),return;end
ids=double(S.(name)(:));
if isfield(S.tie,'ID'), assert(isequal(ids,double(S.tie.ID(:))), 'Tie ID order mismatch'); end
if isfield(S.tie,'id'), assert(isequal(ids,double(S.tie.id(:))), 'Tie ID order mismatch'); end
end

function M=rotation(wpk)
a=deg2rad(wpk(1));b=deg2rad(wpk(2));c=deg2rad(wpk(3));
Rx=[1 0 0;0 cos(a) -sin(a);0 sin(a) cos(a)];
Ry=[cos(b) 0 sin(b);0 1 0;-sin(b) 0 cos(b)];
Rz=[cos(c) -sin(c) 0;sin(c) cos(c) 0;0 0 1];
M=(Rz*Ry*Rx)';
end

function ang=anglesFromRotation(M)
% M' = Rz(kappa)*Ry(phi)*Rx(omega)
A=M';
b=asin(max(-1,min(1,-A(3,1))));
if abs(cos(b))>1e-7
 a=atan2(A(3,2),A(3,3));
 c=atan2(A(2,1),A(1,1));
else
 % Gimbal singularity: choose omega=0 and estimate combined kappa.
 a=0;
 c=atan2(-A(1,2),A(2,2));
end
ang=rad2deg([a b c]);
end

function pix=project(X,pose,IOP)
M=rotation(pose(4:6));
Q=(M*(X-pose(1:3))')';
qz=Q(:,3);
qz(abs(qz)<1e-8)=NaN;
pix=[IOP(1)-IOP(3)*Q(:,1)./qz,IOP(2)+IOP(3)*Q(:,2)./qz];
end

function [pbest,rmse,logs]=solvePose(X,obs,IOP,seeds)
% fminsearch only: no Optimization Toolbox dependency.
% Camera XYZ in mm, rotations in degrees; parameter scaling balances both.
opts=optimset('Display','off','MaxIter',2200,'MaxFunEvals',12000, ...
 'TolX',1e-9,'TolFun',1e-8);
rmse=Inf;pbest=[];logs=zeros(size(seeds,1),3);
for j=1:size(seeds,1)
 p0=seeds(j,:);
 y0=[p0(1:3)/1000,p0(4:6)/30];
 fun=@(v)objective(v,X,obs,IOP);
 [y,fval,~,output]=fminsearch(fun,y0,opts);
 pose=[y(1:3)*1000,y(4:6)*30];
 pred=project(X,pose,IOP);
 d=pred-obs; err=sqrt(mean(sum(d.^2,2)));
 if any(~isfinite(d(:))),err=Inf;end
 logs(j,:)=[j,err,output.iterations];
 if err<rmse,rmse=err;pbest=pose;end
end
end

function v=objective(y,X,obs,IOP)
pose=[y(1:3)*1000,y(4:6)*30];
P=project(X,pose,IOP);
if any(~isfinite(P(:))) || any(abs(P(:))>1e8)
 v=1e20;return;
end
e=P-obs;
% Mild limitation on unphysical camera distance, preventing runaway poses.
v=mean(sum(e.^2,2))+max(0,norm(pose(1:3))-30000)^2/100;
end

function writeRows(file,A,fmt)
f=fopen(file,'w');assert(f>=0,'Cannot write %s',file);
c=onCleanup(@()fclose(f));
for k=1:size(A,1),fprintf(f,fmt,A(k,:));end
end
