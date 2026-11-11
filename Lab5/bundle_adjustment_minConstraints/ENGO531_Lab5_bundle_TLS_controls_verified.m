function ENGO531_Lab5_bundle_TLS_controls_verified()
% ENGO 531 Lab 5: scanner-datum bundle using full_matlab_bundle.p.
% Verified input syntax: ENGO 531 Lab 1 handout, pp. 2-4.
% MODE 0 uses THREE OR MORE *WEIGHTED* control-coordinate observations.
% They are NOT exact fixed constants: sd_XYZ controls their weighting.
% No EOP parameter observations, no measured distances; fixed IOPs.
%
% Inputs: Lab 4 *_results.mat; TLS-initial .ext, .tie; ALL_TLS_controls.csv.
% Rename the supplied P-code to full_matlab_bundle.p and add to path.
% Outputs go in a NEW folder and will not overwrite Lab 4 files.
clc;
FIXED_NAMES = ["F1";"L11";"R12"]; % change controls here
SD_IMAGE_PX = 1.0;
SD_TLS_CONTROL_MM = 2.0;  % trial precision; choose from TLS point-picking quality
SD_UNUSED_DISTANCE_MM = 2.0; % positive, required CFG field but no .dis
SD_UNUSED_EOP_ANGLE_DEG = 0.01; % required CFG field, unused in mode 0

[mf,md] = uigetfile('*.mat','Select FINAL Lab 4 *_results.mat');
if isequal(mf,0),return;end
S=load(fullfile(md,mf));
assert(all(isfield(S,{'names','numericID','photoTable'})), ...
  'Lab 4 MAT file needs names, numericID, photoTable');
names=upper(strtrim(string(S.names(:)))); IDs=double(S.numericID(:));
pho=double(S.photoTable);
assert(numel(names)==numel(IDs) && numel(unique(IDs))==numel(IDs), ...
  'Invalid name/ID mapping');
assert(size(pho,2)==4 && all(isfinite(pho(:))), ...
  'PHO must be [pointID imageID x y]');
assert(all(ismember(pho(:,1),IDs)),'PHO point IDs missing from Lab 4 mapping');
assert(size(unique(pho(:,1:2),'rows'),1)==size(pho,1), ...
  'Duplicate point/image observation');

[ef,ed]=uigetfile('*.ext','Select lab5_tls_initial.ext');
if isequal(ef,0),return;end
ext=readmatrix(fullfile(ed,ef),'FileType','text');
assert(size(ext,2)==8 && all(isfinite(ext(:))), ...
  'EXT must be [imageID cameraID Xc Yc Zc omega phi kappa]');
assert(isequal(sort(ext(:,1)),sort(unique(pho(:,2)))), ...
  'EXT images must exactly match image IDs in PHO');

[tf,td]=uigetfile('*.tie','Select lab5_tls_initial.tie');
if isequal(tf,0),return;end
tieInitial=readmatrix(fullfile(td,tf),'FileType','text');
assert(size(tieInitial,2)==4 && all(isfinite(tieInitial(:))), ...
  'TIE must be [pointID X Y Z]');
assert(isequal(sort(tieInitial(:,1)),sort(IDs)), ...
  'Initial TIE IDs must match the Lab 4 MAT point IDs');

[cf,cd]=uigetfile('*.csv','Select lab5_tls_initial_all_TLS_controls.csv');
if isequal(cf,0),return;end
C=readtable(fullfile(cd,cf),'TextType','string');
required={'point_id','X_mm','Y_mm','Z_mm','BundlePointID'};
assert(all(ismember(required,C.Properties.VariableNames)), ...
  'Control CSV must contain point_id,X_mm,Y_mm,Z_mm,BundlePointID');
cn=upper(strtrim(string(C.point_id)));
assert(numel(unique(cn))==numel(cn),'Duplicate control names');
[found,crow]=ismember(FIXED_NAMES,cn);
assert(all(found),'Some selected controls absent from the TLS CSV');
[found,pRow]=ismember(FIXED_NAMES,names);
assert(all(found),'Some selected controls absent from Lab 4 network');
controlIDs=IDs(pRow);
assert(isequal(controlIDs,double(C.BundlePointID(crow))), ...
  'Control point IDs do not match Lab 4 mapping');
controlXYZ=double([C.X_mm(crow),C.Y_mm(crow),C.Z_mm(crow)]);
assert(all(isfinite(controlXYZ(:))),'Nonfinite scanner control XYZ');
triangleArea=norm(cross(controlXYZ(2,:)-controlXYZ(1,:), ...
                              controlXYZ(3,:)-controlXYZ(1,:)))/2;
assert(triangleArea>1e5,'Selected control triangle is nearly degenerate');

% CRITICAL: one point cannot appear in both .CON and .TIE.
% Controls are represented ONLY by .CON. Uncontrolled points remain in .TIE.
conRows=[controlIDs controlXYZ];
tieRows=tieInitial(~ismember(tieInitial(:,1),controlIDs),:);
assert(isempty(intersect(tieRows(:,1),conRows(:,1))),'Duplicate CON/TIE IDs');
assert(isequal(sort([tieRows(:,1);conRows(:,1)]),sort(IDs)), ...
  'Control/tie partition does not cover all measured object points');

if isfield(S,'IOP')
    IOP=double(S.IOP(:)');
else
    IOP=[2834.0 2156.7 3995.9];
    warning('IOP absent from Lab 4 MAT: using script fallback [xp yp c].');
end
assert(numel(IOP)==3 && all(isfinite(IOP)),'IOP should be [xp yp c]');
intRow=[99 IOP -1]; % fixed, left-handed pixel coordinates
% CFG: mode=0 observed control XYZ, iop_mode=0 fixed IOP.
% All 9 CFG entries must be positive except the two mode flags.
cfg=[0;0;0.001;0.0001;0.01;SD_IMAGE_PX; ...
     SD_UNUSED_DISTANCE_MM;SD_TLS_CONTROL_MM;SD_UNUSED_EOP_ANGLE_DEG];

base=uigetdir(md,'Select output parent folder for Lab 5');
if isequal(base,0),return;end
outDir=fullfile(base,'lab5_TLS_mode0_bundle');
if ~isfolder(outDir),mkdir(outDir);end
prefix=fullfile(outDir,'lab5_TLS');
writeRows([prefix '.pho'],pho,'%d %d %.4f %.4f\n');
writeRows([prefix '.int'],intRow,'%d %.8f %.8f %.8f %d\n');
writeRows([prefix '.ext'],ext,'%d %d %.9f %.9f %.9f %.9f %.9f %.9f\n');
writeRows([prefix '.tie'],tieRows,'%d %.9f %.9f %.9f\n');
writeRows([prefix '.con'],conRows,'%d %.9f %.9f %.9f\n');
writeRows([prefix '.cfg'],cfg,'%g\n');
if isfile([prefix '.dis']),delete([prefix '.dis']);end
writetable(table(names,IDs,'VariableNames',{'point_id','BundlePointID'}), ...
  [prefix '_point_map.csv']);
writetable(table(FIXED_NAMES,controlIDs,controlXYZ(:,1), ...
  controlXYZ(:,2),controlXYZ(:,3), ...
  'VariableNames',{'point_id','BundlePointID','X_mm','Y_mm','Z_mm'}), ...
  [prefix '_control_observations.csv']);

nImages=size(ext,1); nTie=size(tieRows,1); nCon=size(conRows,1);
nPhoto=size(pho,1);
% Mode 0: controls are unknown XYZ, OBSERVED via 3*nCon coordinates.
nUnknown=6*nImages+3*nTie+3*nCon;
nObservations=2*nPhoto+3*nCon;
dof=nObservations-nUnknown;
fprintf('\nLAB 5 MODE 0: %d images, %d ties, %d controls, %d image points\n', ...
  nImages,nTie,nCon,nPhoto);
fprintf('Control names: %s; control triangle area %.3f m^2\n', ...
  strjoin(cellstr(FIXED_NAMES),', '),triangleArea/1e6);
fprintf('n=%d, u=%d, redundancy=%d\n',nObservations,nUnknown,dof);
fprintf('TLS controls are WEIGHTED observations (SD %.3f mm), not hard-fixed.\n', ...
  SD_TLS_CONTROL_MM);
fprintf('Inputs: %s.{pho,int,ext,tie,con,cfg}; NO .dis\n',prefix);
assert(dof>0,'Nonpositive degrees of freedom');

assert(exist('full_matlab_bundle','file')==2 || ...
       exist('full_matlab_bundle','file')==6, ...
       'Rename the provided P-code to full_matlab_bundle.p and add its folder to MATLAB path');
[eop,iop,con,tie,Cx,vr_img,vr_eop,vr_dis,vr_con]=full_matlab_bundle(prefix);
% In mode 0 the course describes vr_eop and vr_dis as ZERO when absent,
% not necessarily MATLAB empty matrices.
assert(noResiduals(vr_eop),'Unexpected EOP observation residuals');
assert(noResiduals(vr_dis),'Unexpected distance residuals');
assert(isstruct(eop) && isfield(eop,'XYZc') && isfield(eop,'wpk'), ...
  'Unexpected EOP output structure');
assert(isstruct(tie) && isfield(tie,'XYZ') && ...
  size(tie.XYZ,1)==nTie,'Unexpected tie output structure');
assert(isstruct(con) && isfield(con,'XYZ') && ...
  size(con.XYZ,1)==nCon,'Unexpected control output structure');
assert(size(vr_img,1)==2*nPhoto,'Image residual count mismatch');
assert(size(vr_con,1)==3*nCon,'Control residual count mismatch');

% Do not assert control XYZ is unchanged: the software treats it as measured
% with sd_XYZ; the estimated control XYZ can legitimately move slightly.
% Lab 1 output order follows EOP, tie, control; verify point ID fields if
% present before interpreting an individual row.
conEstimate=con.XYZ;
conDifference=conEstimate-controlXYZ;
maxControlShift=max(abs(conDifference(:)));
if maxControlShift>3*SD_TLS_CONTROL_MM
    warning('Largest control coordinate correction %.3f mm exceeds 3 * input SD',maxControlShift);
end
imageResiduals=reshape(vr_img(:,1),2,[])';
imageRMS=sqrt(mean(vr_img(:,1).^2));
controlResidualRMS=sqrt(mean(vr_con(:,1).^2));

% Nominal variance factor based on the two ACTUAL observation groups.
vPvImage=sum((vr_img(:,1)/SD_IMAGE_PX).^2);
vPvControl=sum((vr_con(:,1)/SD_TLS_CONTROL_MM).^2);
varianceFactor=(vPvImage+vPvControl)/dof;
fprintf('Image RMS %.4f pixels; control residual RMS %.4f mm\n', ...
 imageRMS,controlResidualRMS);
fprintf('Estimated variance factor %.6f (not automatically applied to Cx)\n',varianceFactor);
fprintf('Max component control adjustment %.4f mm\n',maxControlShift);

writetable(table(pho(:,1),pho(:,2),pho(:,3),pho(:,4), ...
    imageResiduals(:,1),imageResiduals(:,2),vr_img(1:2:end,2), ...
    vr_img(2:2:end,2), ...
    'VariableNames',{'BundlePointID','ImageID','x_pix','y_pix', ...
      'vx_pix','vy_pix','rx','ry'}),[prefix '_image_residuals.csv']);
writetable(table(ext(:,1),eop.XYZc(:,1),eop.XYZc(:,2),eop.XYZc(:,3), ...
    eop.wpk(:,1),eop.wpk(:,2),eop.wpk(:,3), ...
    'VariableNames',{'ImageID','Xc_mm','Yc_mm','Zc_mm','omega_deg', ...
      'phi_deg','kappa_deg'}),[prefix '_adjusted_EOPs.csv']);
writetable(table(tieRows(:,1),tie.XYZ(:,1),tie.XYZ(:,2),tie.XYZ(:,3), ...
    'VariableNames',{'BundlePointID','X_mm','Y_mm','Z_mm'}), ...
    [prefix '_adjusted_ties.csv']);
writetable(table(FIXED_NAMES,controlIDs,controlXYZ(:,1),controlXYZ(:,2), ...
    controlXYZ(:,3),conEstimate(:,1),conEstimate(:,2),conEstimate(:,3), ...
    conDifference(:,1),conDifference(:,2),conDifference(:,3), ...
    'VariableNames',{'point_id','BundlePointID','Xobserved_mm','Yobserved_mm', ...
      'Zobserved_mm','Xestimated_mm','Yestimated_mm','Zestimated_mm', ...
      'dX_mm','dY_mm','dZ_mm'}),[prefix '_control_comparison.csv']);
save([prefix '_results.mat'],'eop','iop','con','tie','Cx', ...
 'vr_img','vr_eop','vr_dis','vr_con','conRows','tieRows','pho','ext', ...
 'cfg','FIXED_NAMES','names','IDs','imageRMS','controlResidualRMS', ...
 'varianceFactor','dof','nUnknown','nObservations');
fprintf('Saved outputs to %s\n',outDir);
fprintf('IMPORTANT: confirm the routine reported CONVERGENCE in the Command Window.\n');
end

function tf=noResiduals(R)
tf=isempty(R) || (isnumeric(R) && isscalar(R) && R==0);
end
function writeRows(path,A,fmt)
fid=fopen(path,'w');assert(fid>=0,'Unable to write %s',path);
c=onCleanup(@()fclose(fid));
for j=1:size(A,1),fprintf(fid,fmt,A(j,:));end
end
