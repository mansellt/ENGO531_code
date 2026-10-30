function ENGO531_Lab4_joint_triangulated()
% ENGO531 Lab 4: prepare ONE 9-image adjustment from the two converged MATs.
% Corrects the key issue in earlier scripts: point IDs 1001... are recycled
% independently in the six-image and three-image bundles. Matching is by
% their saved INITIAL XYZ coordinates against the all-image .tie file.
% Prerequisites: two converged *_results.mat files, all-image .tie file,
% updated calibration(4).PHO file, and full_matlab_bundle.p on MATLAB path.
% Never edits source files. Uses fixed calibrated IOP, Image 1 weighted datum,
% and F1--F6 = 2478 mm as the ONLY scale constraint in the final joint test.
clc;
[f,d]=uigetfile('*.mat','Select CONVERGED images 1-6 results MAT'); if isequal(f,0),return;end
A=load(fullfile(d,f));
[f,d]=uigetfile('*.mat','Select CONVERGED images 7-9 results MAT'); if isequal(f,0),return;end
B=load(fullfile(d,f));
[f,d]=uigetfile('*.tie','Select all-nine original .tie file'); if isequal(f,0),return;end
Tall=readmatrix(fullfile(d,f),'FileType','text');
[f,d]=uigetfile({'*.PHO;*.pho','Raw point observations'},'Select updated calibration(4).PHO'); if isequal(f,0),return;end
fid=fopen(fullfile(d,f),'r'); c=textscan(fid,'%s%f%f%f');fclose(fid);
rawName=string(c{1}); rawIm=c{2}; rawXY=[c{3} c{4}];
folder=uigetdir(d,'Select output folder for JOINT files'); if isequal(folder,0),return;end
prefix=fullfile(folder,'ENGO531_JOINT_1to9');
assert(size(A.eop.XYZc,1)==6 && size(B.eop.XYZc,1)==3,'Select the correct two MAT files');
assert(isfield(A,'tieTable')&&isfield(B,'tieTable')&&isfield(A,'photoTable')&&isfield(B,'phoTable'), ...
 'Missing saved input tables. These exact MAT files should have them.');
assert(size(Tall,2)==4 && size(Tall,1)>=29,'Original all-image .tie seems incomplete');
assert(size(A.tie.XYZ,1)==size(A.tieTable,1) && size(B.tie.XYZ,1)==size(B.tieTable,1), ...
 'Saved tie result differs from saved input rows.');

% Recover identity from actual input 3D approximations, not numericIDs.
% These are shared across the original all-image tie file and each subset's
% saved tieTable; allow 10 mm for R4/R1 length 2701 vs 2704 discrepancy.
ixA=matchXYZ(A.tieTable(:,2:4),Tall(:,2:4),1.0,'images 1-6');
ixB=matchXYZ(B.tieTable(:,2:4),Tall(:,2:4),10,'images 7-9');
fprintf('Matched %d left/front points and %d right-wall points to 9-image initial geometry.\n',numel(ixA),numel(ixB));
assert(isempty(intersect(ixA,ixB)),'Unexpected common point match between distinct faces');
% For the original all-nine input .tie the IDs are assigned to the stable
% first-appearance order of point names in the raw PHO, after removing
% targets with fewer than 2 images. Reconstruct and verify by the X/Y/Z shape.
rawOK=ismember(rawIm,1:9);rawName=rawName(rawOK);rawIm=rawIm(rawOK);rawXY=rawXY(rawOK,:);
[uniqueNames,~,group]=unique(rawName,'stable');count=accumarray(group,1);
retained=uniqueNames(count>=2);
assert(numel(retained)==size(Tall,1),'All-image .tie point count does not match raw PHO filtered to >=2 images.');
expectedID=(1001:1000+numel(retained))';
assert(isequal(Tall(:,1),expectedID),'All-image .tie target IDs are not the expected sequential IDs.');
[isFound,groupFiltered]=ismember(rawName,retained);
P=[expectedID(groupFiltered(isFound)),rawIm(isFound),rawXY(isFound,:)];
assert(all(ismember(1:9,P(:,2))),'A camera has no retained measurements.');

% Validated collinearity convention from geometry test:
% M = (Rx(omega)*Ry(phi)*Rz(kappa))'  (world to camera)
% x = xp - c*qx/qz; y = yp + c*qy/qz; camera y-axis direction -1.
IOP=[2834 2156.7 3995.9];
if isfield(A,'IOP'),IOP=double(A.IOP(:)');end
assert(norm(IOP-double(B.IOP(:)'))<0.1,'The two independent adjustments use different IOPs');
% Explicit validation of the six-image MAT, BEFORE joining the networks.
PA=A.photoTable;[ok,kA]=ismember(PA(:,1),A.tieTable(:,1));assert(all(ok));
err=[];
fprintf('\nREFERENCE REPROJECTION, converged 1-6 solution:\n');
for i=1:6
 mask=PA(:,2)==i;
 predicted=project(A.tie.XYZ(kA(mask),:),[A.eop.XYZc(i,:) A.eop.wpk(i,:)],IOP);
 e=hypot(predicted(:,1)-PA(mask,3),predicted(:,2)-PA(mask,4));
 fprintf(' Image %d: %.3f px radial RMS (%d points)\n',i,sqrt(mean(e.^2)),numel(e));
 if i>=2,err=[err;e];end %#ok<AGROW>
end
refRMS=sqrt(mean(err.^2));
fprintf(' Verified images 2-6 reference RMS = %.3f px\n',refRMS);
assert(isfinite(refRMS)&&refRMS<15,'Reference validation failed; do not run joint bundle.');

% Resection: solve camera 6 in the right-wall coordinate frame from R targets
% observed in original Image 6 (not present in the 7-9 MAT's photoTable).
[ok,loc]=ismember(Tall(ixB,1),Tall(:,1));assert(all(ok));
[has,j]=ismember(P(:,1),Tall(ixB,1));
sel=P(:,2)==6 & has;
obs6=P(sel,3:4);idxB=j(sel);
assert(numel(idxB)>=8,'Cannot find enough image-6 R-wall observations');
Xr=B.tie.XYZ(idxB,:);
seeds=[3486 -1652 -495 90 30 0;3942 -933 170 67 45 26; ...
       4000 -1300 -300 90 45 0;4200 -500 -500 90 55 0; ...
       B.eop.XYZc(1,:) B.eop.wpk(1,:)];
opts=optimset('Display','off','MaxIter',4000,'MaxFunEvals',25000,'TolX',1e-8,'TolFun',1e-8);
best=Inf;bestP=[];
for ii=1:size(seeds,1)
 init=[seeds(ii,1:3)/1000,seeds(ii,4:6)/30];
 [v,~,flag]=fminsearch(@(v)poseCost(v,Xr,obs6,IOP),init,opts);
 p=[v(1:3)*1000 v(4:6)*30];
 pr=project(Xr,p,IOP);
 er=sqrt(mean(sum((pr-obs6).^2,2)));
 fprintf('Image 6 resection seed %d: RMS %.3f px; exit flag %d\n',ii,er,flag);
 if isfinite(er) && er<best,best=er;bestP=p;end
end
fprintf('Best Image 6 resection RMS: %.3f px\n',best);
if isempty(bestP)||best>35
 fprintf(2,'STOP: Image 6 resection not accurate enough to register networks. No bundle run.\n');
 return;
end

% Frame transform: worldLeft = R*worldRight + translation.
ML=rotW2C(A.eop.wpk(6,:)); MR=rotW2C(bestP(4:6));
R=ML'*MR; translation=A.eop.XYZc(6,:)'-R*bestP(1:3)';
TallNew=Tall;
% For points missing from the saved 1-6 solution, intersect camera rays
% from the ORIGINAL image observations using the CONVERGED 1-6 EOPs.
% This is essential: a rigid transform of estimated face coordinates is
% not precise enough for the extra L/F points in Images 1-3.
extraLF=setdiff((1:size(Tall,1))',[ixA;ixB]);
assert(all(startsWith(retained(extraLF),"L")|startsWith(retained(extraLF),"F")), ...
    'An unexpected unmatched point is not an L/F target.');
fprintf('\nTRIANGULATING %d missing L/F points from calibrated camera rays:\n',numel(extraLF));
triResults=zeros(numel(extraLF),4);
for j=1:numel(extraLF)
    idx=extraLF(j); pid=Tall(idx,1);
    mask=P(:,1)==pid & P(:,2)<=6;
    assert(nnz(mask)>=2,'Target %s has fewer than two 1-6 observations',retained(idx));
    obs=P(mask,:);
    [X,linearRank]=triangulateWeightedLinear(obs,A.eop.XYZc,A.eop.wpk,IOP);
    assert(linearRank==3 && all(isfinite(X)), ...
        'Triangulation failed for target %s. Inspect stereo geometry.',retained(idx));
    % Refine 3D point alone, with all six-image EOPs FIXED.
    options=optimset('Display','off','MaxIter',1500,'MaxFunEvals',5000, ...
        'TolX',1e-9,'TolFun',1e-9);
    X=fminsearch(@(x) pointCost(x,obs,A.eop.XYZc,A.eop.wpk,IOP),X,options);
    errs=pointErrors(X,obs,A.eop.XYZc,A.eop.wpk,IOP);
    rms=sqrt(mean(errs.^2));
    fprintf('  %-5s images %s: radial RMS %8.3f px, XYZ [%9.2f %9.2f %9.2f] mm\n', ...
        char(retained(idx)),mat2str(obs(:,2)'),rms,X(1),X(2),X(3));
    TallNew(idx,2:4)=X;
    triResults(j,:)=[pid,nnz(mask),rms,max(errs)];
end
TallNew(ixA,2:4)=A.tie.XYZ;
TallNew(ixB,2:4)=(R*B.tie.XYZ'+translation)';
if any(triResults(:,3)>15)
    warning('Some new points have >15px triangulation RMS. Inspect corresponding image labels.');
end
writetable(array2table(triResults,'VariableNames', ...
  {'ID','NumberOfRays','RadialRMS_px','MaxError_px'}),[prefix '_triangulation.csv']);
% Original source EOPs for image numbering; replace with converged poses.
notes=[1 -2213 1711 -645 -90 0 90;2 -2279 561 -495 -75 0 90; ...
3 -1619 -1540 -495 -45 0 90;4 -440 -2027 -495 -25 0 90; ...
5 1239 -2545 -495 0 0 90;6 3486 -1652 -495 30 0 90; ...
7 4368 -97 -495 60 0 90;8 4901 411 -495 90 0 90;9 4455 2686 -495 120 0 90];
E=[notes(:,1),99*ones(9,1),notes(:,2:4),notes(:,7),notes(:,5),notes(:,6)];
E(1:6,3:5)=A.eop.XYZc; E(1:6,6:8)=A.eop.wpk;
for i=7:9
 k=i-6;
 E(i,3:5)=(R*B.eop.XYZc(k,:)'+translation)';
 Mnew=rotW2C(B.eop.wpk(k,:))*R';
 E(i,6:8)=rot2WPK(Mnew);
end
% Preflight per camera against ORIGINAL measurements (raw PHO not changed)
[ok,k]=ismember(P(:,1),TallNew(:,1));assert(all(ok));
summary=zeros(9,3);
for i=1:9
 mask=P(:,2)==i; pred=project(TallNew(k(mask),2:4),E(i,3:8),IOP);
 e=hypot(pred(:,1)-P(mask,3),pred(:,2)-P(mask,4));
 summary(i,:)=[i,nnz(mask),sqrt(mean(e.^2))];
end
fprintf('\nPRE-BUNDLE REPROJECTION RMS, each image:\n');
disp(array2table(summary,'VariableNames',{'Image','Count','RadialRMS_px'}));
writetable(array2table(summary,'VariableNames',{'Image','Count','RadialRMS_px'}),[prefix '_preflight.csv']);
if max(summary(:,3))>100 || any(~isfinite(summary(:,3)))
 fprintf(2,'STOP: initial geometry is not consistent enough. Joint bundle not run.\n');
 return;
end
% Create NEW course files, original image observations preserved.
scaleA=find(retained=="F1",1);scaleB=find(retained=="F6",1);
assert(~isempty(scaleA)&&~isempty(scaleB),'Scale endpoints F1/F6 absent');
writeRows([prefix '.pho'],P,'%d %d %.4f %.4f\n');
writeRows([prefix '.tie'],TallNew,'%d %.7f %.7f %.7f\n');
writeRows([prefix '.ext'],E,'%d %d %.7f %.7f %.7f %.8f %.8f %.8f\n');
writeRows([prefix '.int'],[99 IOP -1],'%d %.5f %.5f %.5f %d\n');
r4 = find(retained == "R4");
r1 = find(retained == "R1");

distanceObservations = [
    expectedID(scaleA), expectedID(scaleB), 2478;
    expectedID(r4),     expectedID(r1),     2701
];

writeRows([prefix '.dis'], distanceObservations, ...
    '%d %d %.4f\n');
cfg = [1;0;0.01;0.001;0.01;3.2;5;1;0.01];
writeRows([prefix '.cfg'],cfg,'%g\n');
writetable(table(retained,expectedID,'VariableNames',{'Target','ID'}),[prefix '_point_map.csv']);
fprintf('Prepared ONE nine-image adjustment: %s.*\n',prefix);
reply=input('Run the nine-image course adjustment? [y/N]: ','s');
if ~strcmpi(strtrim(reply),'y'),return;end
[eop,iop,con,tie,Cx,vr_img,vr_eop,vr_dis,vr_con]=full_matlab_bundle(prefix);
fprintf('Returned RMS: %.3f px (valid ONLY if course explicitly reports convergence)\n',sqrt(mean(vr_img(:,1).^2)));
save([prefix '_results.mat'],'eop','iop','con','tie','Cx','vr_img','vr_eop','vr_dis','vr_con','E','TallNew','R','translation','summary');
end
function [X,rk]=triangulateWeightedLinear(obs,cameras,angles,IOP)
% Each observation gives two linearized collinearity equations:
% (m1 - xnorm*m3)*(X-C)=0; (m2 - ynorm*m3)*(X-C)=0
n=size(obs,1); AA=zeros(2*n,3); bb=zeros(2*n,1);
for k=1:n
    im=obs(k,2); M=rotW2C(angles(im,:)); C=cameras(im,:)';
    xn=-(obs(k,3)-IOP(1))/IOP(3);
    yn= (obs(k,4)-IOP(2))/IOP(3);
    G=[M(1,:)-xn*M(3,:);M(2,:)-yn*M(3,:)];
    AA(2*k-1:2*k,:)=G; bb(2*k-1:2*k)=G*C;
end
rk=rank(AA);X=(AA\bb)';
end
function errors=pointErrors(X,obs,cameras,angles,IOP)
errors=zeros(size(obs,1),1);
for k=1:size(obs,1)
  im=obs(k,2); proj=project(X,[cameras(im,:) angles(im,:)],IOP);
  errors(k)=norm(proj-obs(k,3:4));
end
end
function cost=pointCost(X,obs,cameras,angles,IOP)
e=pointErrors(X,obs,cameras,angles,IOP);
if any(~isfinite(e)) || max(abs(e))>1e7
 cost=1e18;
else
 cost=mean(e.^2);
end
end
function [R,t,rms]=rigidFit(initial,adjusted)
center0=mean(initial,1)';center1=mean(adjusted,1)';
H=(initial-center0')'*(adjusted-center1');
[U,~,V]=svd(H);
R=V*diag([1,1,det(V*U')])*U';
t=center1-R*center0;
d=(R*initial'+t)-adjusted';
rms=sqrt(mean(sum(d.^2,1)));
end
function ix=matchXYZ(src,dst,maxError,tag)
D=zeros(size(src,1),size(dst,1));
for j=1:size(src,1),D(j,:)=sqrt(sum((dst-src(j,:)).^2,2));end
[val,ix]=min(D,[],2);
assert(all(val<maxError),'Some %s tie coordinates have no matching all-nine initialization',tag);
assert(numel(unique(ix))==numel(ix),'Ambiguous %s tie matching',tag);
end
function M=rotW2C(wpk)
a=deg2rad(wpk(1)); b=deg2rad(wpk(2)); c=deg2rad(wpk(3));
Rx=[1 0 0;0 cos(a) -sin(a);0 sin(a) cos(a)];
Ry=[cos(b) 0 sin(b);0 1 0;-sin(b) 0 cos(b)];
Rz=[cos(c) -sin(c) 0;sin(c) cos(c) 0;0 0 1];
M=(Rx*Ry*Rz)';
end
function ang=rot2WPK(M)
A=M';b=asin(max(-1,min(1,A(1,3))));
if abs(cos(b))>1e-7
 a=atan2(-A(2,3),A(3,3));c=atan2(-A(1,2),A(1,1));
else
 a=atan2(A(3,2),A(2,2)); c=0;
end
ang=rad2deg([a b c]);
end
function pix=project(X,p,IOP)
M=rotW2C(p(4:6));Q=(M*(X-p(1:3))')';
z=Q(:,3);z(abs(z)<1e-9)=NaN;
pix=[IOP(1)-IOP(3)*Q(:,1)./z,IOP(2)+IOP(3)*Q(:,2)./z];
end
function v=poseCost(x,X,obs,IOP)
p=[x(1:3)*1000 x(4:6)*30];pred=project(X,p,IOP);
if any(~isfinite(pred(:)))||any(abs(pred(:))>1e8),v=1e20;return;end
d=pred-obs;v=mean(sum(d.^2,2))+max(0,norm(p(1:3))-20000)^2/200;
end
function writeRows(file,A,fmt)
f=fopen(file,'w');assert(f>=0,'Cannot write output');c=onCleanup(@()fclose(f));
for k=1:size(A,1),fprintf(f,fmt,A(k,:));end
end
