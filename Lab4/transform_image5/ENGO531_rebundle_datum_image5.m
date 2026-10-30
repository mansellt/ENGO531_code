function ENGO531_rebundle_datum_image5()
% Re-run the ONE nine-photo bundle with ORIGINAL IMAGE 5 as the datum.
% Uses the converged nine-photo output as starting geometry, not field notes.
% Needs these files in one folder:
% ENGO531_JOINT_1to9_results.mat, .pho, .int, .dis, .cfg
% and full_matlab_bundle.p on your MATLAB path.
% Uses the same point observations, two tape distances, and fixed IOPs.
clc;
[fn,fd] = uigetfile('*.mat','Choose the converged ENGO531_JOINT_1to9_results.mat');
if isequal(fn,0),return;end
S=load(fullfile(fd,fn));
assert(isfield(S,'eop') && isfield(S,'tie'), 'MAT must contain adjusted eop and tie.');
assert(size(S.eop.XYZc,1)==9 && size(S.eop.wpk,1)==9, 'Expected nine adjusted EOPs.');
assert(isfield(S,'TallNew'), 'Expected TallNew point-ID mapping from standalone joint script.');
assert(size(S.tie.XYZ,1)==size(S.TallNew,1), 'Adjusted tie rows do not match TallNew.');
[~,base,~]=fileparts(fn);
base=regexprep(base,'_results$','');
phoFile=fullfile(fd,[base '.pho']); disFile=fullfile(fd,[base '.dis']);
intFile=fullfile(fd,[base '.int']); cfgFile=fullfile(fd,[base '.cfg']);
assert(isfile(phoFile) && isfile(disFile) && isfile(intFile) && isfile(cfgFile), ...
 'Place the original .pho/.dis/.int/.cfg files next to the results MAT.');
pho=readmatrix(phoFile,'FileType','text');
dis=readmatrix(disFile,'FileType','text');
int=readmatrix(intFile,'FileType','text');
cfg=readmatrix(cfgFile,'FileType','text');
pho=pho(:,1:4); dis=dis(:,1:3); int=int(1,1:5); cfg=cfg(:);
assert(all(ismember(unique(pho(:,2)),1:9)), 'Source photo image IDs differ from 1:9.');
assert(numel(unique(pho(:,2)))==9,'Missing image observations.');
% Point rows in S.tie.XYZ follow S.TallNew, verified in standalone setup.
tie=[double(S.TallNew(:,1)),double(S.tie.XYZ)];
C=double(S.eop.XYZc); angles=double(S.eop.wpk);
% Set true to express XYZ and camera EOPs in the axes of Image 5.
% The camera-5 axes are NOT necessarily the object's face axes.
alignAxesWithImage5=true;
if alignAxesWithImage5
    M5=rotW2C(angles(5,:)); O=C(5,:)';
    xyz=(M5*(tie(:,2:4)'-O))';
    CC=(M5*(C'-O))';
    AA=zeros(9,3);
    for i=1:9
       Mnew=rotW2C(angles(i,:))*M5';
       AA(i,:)=rot2WPK(Mnew);
    end
    % Force exact identity for datum image 5 (numerically equivalent).
    CC(5,:)=[0 0 0]; AA(5,:)=[0 0 0];
else
    xyz=tie(:,2:4);CC=C;AA=angles;
end
% Validate that re-datuming changes no predicted image coordinates.
IOP=int(2:4);[ok,loc]=ismember(pho(:,1),tie(:,1));
assert(all(ok),'An observation points to a missing target ID.');
maxPixDiff=0;
for i=1:9
  rows=pho(:,2)==i;
  pOld=project(C(i,:),angles(i,:),tie(loc(rows),2:4),IOP);
  pNew=project(CC(i,:),AA(i,:),xyz(loc(rows),:),IOP);
  dif=hypot(pOld(:,1)-pNew(:,1),pOld(:,2)-pNew(:,2));
  maxPixDiff=max(maxPixDiff,max(dif));
end
fprintf('Maximum change in predicted image coordinates = %.6g px\n',maxPixDiff);
assert(maxPixDiff<0.01,'Re-datum conversion changed the projection. Stop.');
% The course constrains the FIRST CAMERA; put original Image 5 first.
origOrder=[5 1 2 3 4 6 7 8 9];
[~,newID]=ismember(pho(:,2),origOrder);
assert(all(newID>0),'Image renumbering failed.');
pho(:,2)=newID;
E=[(1:9)' 99*ones(9,1) CC(origOrder,:) AA(origOrder,:)];
T=[tie(:,1) xyz];
outDir=uigetdir(fd,'Choose a NEW folder for the Image-5 datum run');
if isequal(outDir,0),return;end
prefix=fullfile(outDir,'ENGO531_joint_datum_image5');
assert(~strcmp(outDir,fd) || ~strcmp(base,'ENGO531_joint_datum_image5'), ...
 'Choose a different output name or folder.');
writeRows([prefix '.pho'],pho,'%d %d %.5f %.5f\n');
writeRows([prefix '.tie'],T,'%d %.9f %.9f %.9f\n');
writeRows([prefix '.ext'],E,'%d %d %.9f %.9f %.9f %.10f %.10f %.10f\n');
writeRows([prefix '.int'],int,'%d %.8f %.8f %.8f %d\n');
writeRows([prefix '.dis'],dis,'%d %d %.6f\n');
writeRows([prefix '.cfg'],cfg,'%g\n');
if isfile([prefix '.con']),delete([prefix '.con']);end
writetable(table((1:9)',origOrder','VariableNames',{'BundleImageID','OriginalImageID'}), ...
 [prefix '_image_number_map.csv']);
fprintf('Original Image 5 is bundle Image 1. Image 5 datum EOP: [');
fprintf(' %.5f',E(1,3:8));fprintf(' ]\n');
fprintf('Full network: %d points, %d image observations, %d distances.\n', ...
 size(T,1),size(pho,1),size(dis,1));
fprintf('One bundle input set created: %s.*\n',prefix);
reply=input('Run course adjustment? [y/N]: ','s');
if ~strcmpi(strtrim(reply),'y'),return;end
assert(exist('full_matlab_bundle','file')==2 || exist('full_matlab_bundle','file')==6, ...
 'Add full_matlab_bundle.p to the MATLAB path.');
[eop,iop,con,tieResult,Cx,vr_img,vr_eop,vr_dis,vr_con]=full_matlab_bundle(prefix);
rmsXY=sqrt(mean(vr_img(:,1).^2));
fprintf('Returned RMS = %.4f px (use only if course reports convergence)\n',rmsXY);
printEOP=table(origOrder',eop.XYZc(:,1),eop.XYZc(:,2),eop.XYZc(:,3), ...
 eop.wpk(:,1),eop.wpk(:,2),eop.wpk(:,3), ...
 'VariableNames',{'OriginalImage','Xc_mm','Yc_mm','Zc_mm','Omega_deg','Phi_deg','Kappa_deg'});
disp(printEOP);
save([prefix '_results.mat'],'eop','iop','con','tieResult','Cx', ...
 'vr_img','vr_eop','vr_dis','vr_con','rmsXY','origOrder','alignAxesWithImage5', ...
 'pho','T','E','dis','cfg','int');
end

function M=rotW2C(wpk)
a=deg2rad(wpk(1));b=deg2rad(wpk(2));c=deg2rad(wpk(3));
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
 a=atan2(A(3,2),A(2,2));c=0;
end
ang=rad2deg([a b c]);
end
function pix=project(C,wpk,X,IOP)
M=rotW2C(wpk);Q=(M*(X-C)')';
z=Q(:,3);z(abs(z)<1e-10)=NaN;
pix=[IOP(1)-IOP(3)*Q(:,1)./z,IOP(2)+IOP(3)*Q(:,2)./z];
end
function writeRows(file,A,fmt)
f=fopen(file,'w');assert(f>=0,'Cannot create file: %s',file);
c=onCleanup(@()fclose(f));
for k=1:size(A,1),fprintf(f,fmt,A(k,:));end
end
