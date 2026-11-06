function ENGO531_Lab5_bundle_TLS_controls()
% ENGO 531 LAB 5: Three TLS-fixed points, no distance/EOP observations.
% Calls the SAME course full_matlab_bundle.p as the Lab 4 script.
%
% IMPORTANT INTERFACE CHECK:
% The uploaded .p function does not reveal the .con layout or cfg modes.
% This script ASSUMES .con rows = pointID X Y Z (mm) and cfg(1)=0 means
% no EOP pseudo-observations. Confirm these with the instructor/example
% before setting COURSE_FORMAT_CONFIRMED=true. Until then, files are
% generated but the compiled adjustment is NOT called.
%
% Required inputs:
%   - final Lab 4 *_results.mat (photoTable, names, numericID, etc.)
%   - lab5_tls_initial.ext and lab5_tls_initial.tie (TLS-transformed initials)
%   - lab5_tls_initial_all_TLS_controls.csv (scanner XYZ in millimetres)
%   - full_matlab_bundle.p (can be named full_matlab_bundle(4).p on upload;
%     rename to full_matlab_bundle.p when placing alongside script)
%
% Selection: F1, L11, R12. Change FIXED_NAMES if desired.
% Results written to a new subfolder, without altering Lab 4 files.

clc;
FIXED_NAMES = ["F1";"L11";"R12"];
SD_IMAGE_PX = 1.0;
COURSE_FORMAT_CONFIRMED = false;  % set TRUE only after verifying .con / .cfg

[fn,fd]=uigetfile('*.mat','Select final Lab 4 *_results.mat');
if isequal(fn,0),return;end
S=load(fullfile(fd,fn));
assert(all(isfield(S,{'names','numericID','photoTable'})), ...
    'Lab 4 MAT needs names, numericID, and photoTable.');
names=upper(strtrim(string(S.names(:)))); ids=double(S.numericID(:));
photo=double(S.photoTable);
assert(numel(names)==numel(ids),'names/numericID length mismatch');
assert(size(photo,2)==4,'photoTable must be [pointID,imageID,x_pixel,y_pixel]');
assert(numel(unique(ids))==numel(ids),'Duplicate bundle point IDs');
assert(all(ismember(photo(:,1),ids)),'PHO contains IDs outside the selected Lab 4 network');
assert(numel(unique(photo(:,1)+photo(:,2)/1e5))==size(photo,1), ...
    'Duplicate point/image measurements');

[eName,eDir]=uigetfile('*.ext','Select lab5_tls_initial.ext');
if isequal(eName,0),return;end
ext=readmatrix(fullfile(eDir,eName),'FileType','text');
assert(size(ext,2)==8,'EXT must be [image camera Xc Yc Zc omega phi kappa]');
assert(all(isfinite(ext(:))),'Nonfinite EOP values');
assert(all(ismember(unique(photo(:,2)),ext(:,1))), ...
    'EXT missing at least one active image');
assert(all(ismember(ext(:,1),unique(photo(:,2)))), ...
    'EXT has image IDs not used in the observations');

[tName,tDir]=uigetfile('*.tie','Select lab5_tls_initial.tie');
if isequal(tName,0),return;end
tie=readmatrix(fullfile(tDir,tName),'FileType','text');
assert(size(tie,2)==4,'TIE must be [pointID X Y Z]');
assert(all(isfinite(tie(:))),'Nonfinite tie coordinates');
assert(isequal(sort(tie(:,1)),sort(ids)), ...
    'TIE point IDs do not match MAT numericID list');

[cName,cDir]=uigetfile('*.csv','Select lab5_tls_initial_all_TLS_controls.csv');
if isequal(cName,0),return;end
C=readtable(fullfile(cDir,cName),'TextType','string');
required={'point_id','X_mm','Y_mm','Z_mm','BundlePointID'};
assert(all(ismember(required,C.Properties.VariableNames)), ...
    'TLS controls CSV missing expected columns');
cn=upper(strtrim(string(C.point_id)));
assert(numel(unique(cn))==numel(cn),'Duplicate TLS control names');
[ok,cix]=ismember(FIXED_NAMES,cn);
assert(all(ok),'A requested fixed control is absent from scanner controls');
controlIDs=double(C.BundlePointID(cix));
[ok,idRows]=ismember(FIXED_NAMES,names);
assert(all(ok) && isequal(ids(idRows),controlIDs), ...
    'Control name-to-bundle-ID mapping does not agree with Lab 4 MAT');
controlXYZ=double([C.X_mm(cix),C.Y_mm(cix),C.Z_mm(cix)]);
assert(all(isfinite(controlXYZ(:))),'Nonfinite scanner control coordinates');
area=norm(cross(controlXYZ(2,:)-controlXYZ(1,:), ...
                controlXYZ(3,:)-controlXYZ(1,:)))/2;
assert(area>1e5,'Three controls are too close to collinear');

% Fix the initial XYZ of the chosen control points to their exact scanner
% coordinates, so the solver does not begin with a mismatch at fixed points.
for k=1:3
    row=find(tie(:,1)==controlIDs(k),1);
    tie(row,2:4)=controlXYZ(k,:);
end

% Keep calibrated image parameters fixed. Previous Lab 4 values:
%   99 xp yp c y_dir   (pixels; y_dir=-1)
if isfield(S,'IOP')
    iop=double(S.IOP(:)');
else
    iop=[2834.0,2156.7,3995.9];
    warning('No IOP in MAT: using Lab 4 script constants. Verify Lab 3 calibration.');
end
assert(numel(iop)==3 && all(isfinite(iop)),'Expected fixed [xp yp c]');
int=[99,iop,-1];

% Lab 4: cfg = [1;0;0.001;0.0001;0.01;sd_xy;sd_dist;sd_eop_xyz;sd_eop_angle]
% Here we ASSUME the first entry 0 disables EOP pseudo-observations.
% Keep all nine numeric entries required by the compiled course function.
cfg=[0;0;0.001;0.0001;0.01;SD_IMAGE_PX;2;1;0.01];

base=uigetdir(fd,'Choose output parent folder for Lab 5 bundle');
if isequal(base,0),return;end
outDir=fullfile(base,'lab5_three_tls_fixed_bundle');
if ~isfolder(outDir),mkdir(outDir);end
prefix=fullfile(outDir,'lab5_TLS');

writeNumeric([prefix '.pho'],photo,'%d %d %.4f %.4f\n');
writeNumeric([prefix '.ext'],ext,'%d %d %.9f %.9f %.9f %.9f %.9f %.9f\n');
writeNumeric([prefix '.tie'],tie,'%d %.9f %.9f %.9f\n');
writeNumeric([prefix '.int'],int,'%d %.8f %.8f %.8f %d\n');
writeNumeric([prefix '.cfg'],cfg,'%g\n');
% HYPOTHESIZED course CON format. See note at top; verify BEFORE RUNNING.
controlRows=[controlIDs,controlXYZ];
writeNumeric([prefix '.con'],controlRows,'%d %.9f %.9f %.9f\n');
% No distance file is generated. Delete a stale one if this exact run exists.
if isfile([prefix '.dis']),delete([prefix '.dis']);end
map=table(names,ids,'VariableNames',{'point_id','BundlePointID'});
writetable(map,[prefix '_point_map.csv']);
writetable(table(FIXED_NAMES,controlIDs,controlXYZ(:,1),controlXYZ(:,2), ...
    controlXYZ(:,3),'VariableNames', ...
    {'point_id','BundlePointID','X_mm','Y_mm','Z_mm'}),[prefix '_fixed_points.csv']);

nC=size(ext,1); nP=size(tie,1); nObs=size(photo,1);
% 9 scalar control coordinate conditions (3 control points * 3 dimensions),
% 6 EOP parameters per camera, 3 tie coords per object point.
nominalDOF=2*nObs+9-(6*nC+3*nP);
fprintf('\nLAB 5 input summary:\n');
fprintf('Images: %d; object points: %d; image point measurements: %d\n',nC,nP,nObs);
fprintf('Fixed control triangle area: %.3f m^2\n',area/1e6);
fprintf('Approx. redundancy (9 fixed scalar constraints assumed): %d\n',nominalDOF);
fprintf('Controls: %s\n',strjoin(cellstr(FIXED_NAMES),', '));
fprintf('Written: %s.{pho,ext,tie,int,cfg,con}\n',prefix);
fprintf('No distance file. No first-camera EOP observations requested.\n');
assert(nominalDOF>0,'Network has non-positive nominal redundancy');

if ~COURSE_FORMAT_CONFIRMED
    fprintf(['\nSTOPPED SAFELY BEFORE BUNDLE CALL.\n' ...
      'The uploaded full_matlab_bundle.p is compiled, so .con syntax and\n' ...
      'cfg(1)=0 cannot be established from the supplied source. Verify both\n' ...
      'against course instructions/example. Set COURSE_FORMAT_CONFIRMED=true\n' ...
      'only when verified; then rerun this script.\n']);
    return
end

assert(exist('full_matlab_bundle','file')==2 || ...
       exist('full_matlab_bundle','file')==6, ...
    'Place full_matlab_bundle.p on the MATLAB path (exact filename).');
[eop,iopOut,con,tieOut,Cx,vr_img,vr_eop,vr_dis,vr_con]= ...
    full_matlab_bundle(prefix);

% Reject any result inconsistent with the intended observation model.
assert(isempty(vr_eop),'EOP observations unexpectedly active. CHECK CFG!');
assert(isempty(vr_dis),'Distance observations unexpectedly active. CHECK DIS!');
assert(isfield(eop,'XYZc') && isfield(eop,'wpk'), ...
    'Bundle returned unexpected EOP structure');
assert(isfield(tieOut,'XYZ'),'Bundle returned unexpected tie structure');
assert(size(tieOut.XYZ,1)==nP,'Bundle output point count mismatch');

% Verify fixed coordinates by the same tie-file row order. If the course
% routine orders tieOut differently, DO NOT accept this check until mapped.
if isfield(tieOut,'ID')
    [found,rows]=ismember(controlIDs,double(tieOut.ID(:)));
elseif isfield(tieOut,'id')
    [found,rows]=ismember(controlIDs,double(tieOut.id(:)));
else
    found=true(3,1);
    [~,rows]=ismember(controlIDs,tie(:,1));
    warning(['tieOut lacks explicit IDs: assuming input TIE order is preserved. ' ...
             'Verify before interpreting results.']);
end
assert(all(found),'Cannot match controls in adjusted tie points');
fixedDiff=tieOut.XYZ(rows,:)-controlXYZ;
maxFixedError=max(abs(fixedDiff(:)));
fprintf('Largest fixed-point coordinate discrepancy: %.8g mm\n',maxFixedError);
assert(maxFixedError<1e-4, ...
    'FIXED CONTROLS MOVED: .con constraints not applied as intended!');

v=vr_img(:,1);
assert(numel(v)==2*nObs,'Unexpected number of image residuals');
imageRMS=sqrt(mean(v.^2));
fprintf('Image coordinate RMS: %.5f pixels\n',imageRMS);
if size(Cx,1)>0
    fprintf('Returned covariance/cofactor matrix: %d-by-%d\n',size(Cx,1),size(Cx,2));
end
imageResiduals=reshape(v,2,[])';
residuals=table(photo(:,1),photo(:,2),photo(:,3),photo(:,4), ...
    imageResiduals(:,1),imageResiduals(:,2), ...
    'VariableNames',{'BundlePointID','ImageID','x_pix','y_pix','vx_pix','vy_pix'});
writetable(residuals,[prefix '_image_residuals.csv']);
images=ext(:,1);
eopTable=table(images,eop.XYZc(:,1),eop.XYZc(:,2),eop.XYZc(:,3), ...
    eop.wpk(:,1),eop.wpk(:,2),eop.wpk(:,3), ...
    'VariableNames',{'ImageID','Xc_mm','Yc_mm','Zc_mm','omega_deg','phi_deg','kappa_deg'});
writetable(eopTable,[prefix '_adjusted_EOPs.csv']);
ptTable=table(tie(:,1),tieOut.XYZ(:,1),tieOut.XYZ(:,2),tieOut.XYZ(:,3), ...
    'VariableNames',{'BundlePointID','X_mm','Y_mm','Z_mm'});
writetable(ptTable,[prefix '_adjusted_points.csv']);
save([prefix '_results.mat'],'eop','iopOut','con','tieOut','Cx', ...
    'vr_img','vr_eop','vr_dis','vr_con','FIXED_NAMES','controlIDs', ...
    'controlXYZ','photo','ext','tie','cfg','int','imageRMS', ...
    'fixedDiff','maxFixedError');
fprintf('Saved verified output to %s\n',outDir);
end

function writeNumeric(file,A,fmt)
fid=fopen(file,'w');assert(fid>=0,'Cannot create %s',file);
clean=onCleanup(@()fclose(fid));
for k=1:size(A,1),fprintf(fid,fmt,A(k,:));end
end
