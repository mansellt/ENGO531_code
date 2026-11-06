function prepare_lab5_tls_bundle()
% ENGO 531 Lab 5 - transform Lab 4 approximate network into TLS space.
% Requires the FINAL Lab 4 *_results.mat file and TLS target_coordinates.csv.
% Produces transformed initial .ext and .tie and control-coordinate CSV.
% IMPORTANT: This prepares initial values, NOT the final Lab 5 bundle.
% Your course P-code control-file syntax / CFG modes must be confirmed before
% constructing .con and calling full_matlab_bundle.
%
clc;
[fn,fd] = uigetfile('*.mat','Select FINAL Lab 4 results MAT file');
if isequal(fn,0),return;end
S = load(fullfile(fd,fn));
assert(isfield(S,'names') && isfield(S,'tie') && isfield(S,'eop'), ...
    'Results MAT must contain names, tie, and eop.');
assert(isfield(S.tie,'XYZ') && isfield(S.eop,'XYZc') && isfield(S.eop,'wpk'), ...
    'Expected fields tie.XYZ, eop.XYZc and eop.wpk.');
names = upper(strtrim(string(S.names(:))));
oldTie = double(S.tie.XYZ);
oldCam = double(S.eop.XYZc);
oldWpk = double(S.eop.wpk); % assumed omega,phi,kappa, DEGREES
assert(size(oldTie,1)==numel(names) && size(oldTie,2)==3,'Tie XYZ/name mismatch.');
assert(size(oldCam,2)==3 && isequal(size(oldWpk),size(oldCam)), 'EOP array size mismatch.');

[cn,cd] = uigetfile('*.csv','Select target_coordinates.csv');
if isequal(cn,0),return;end
C = readtable(fullfile(cd,cn),'TextType','string');
assert(all(ismember({'point_id','X','Y','Z'},C.Properties.VariableNames)), ...
    'CSV needs point_id,X,Y,Z.');
controlNames = upper(strtrim(string(C.point_id)));
assert(numel(unique(controlNames)) == numel(controlNames),'Duplicate TLS target IDs.');
[matched,loc] = ismember(controlNames,names);
assert(nnz(matched)>=3, 'Need at least 3 matching TLS and Lab 4 points.');
oldMatch = oldTie(loc(matched),:);
newMatch = 1000*double([C.X(matched),C.Y(matched),C.Z(matched)]); % TLS m -> mm
assert(all(isfinite(oldMatch(:))) && all(isfinite(newMatch(:))), 'Nonfinite control coordinates.');

% Similarity fit: new(mm) = scale * old(mm) * Rrow + translation.
% Kabsch SVD enforces a proper rotation (not a reflection).
a = mean(oldMatch,1); b = mean(newMatch,1);
A = oldMatch-a; B = newMatch-b;
assert(rank(A,1e-7)>=2,'Matched control points are collinear.');
[U,~,V] = svd(A'*B);
D=diag([1 1 sign(det(U*V'))]);
Rrow=U*D*V';
scale=sum(sum((A*Rrow).*B))/sum(A(:).^2);
assert(isfinite(scale) && scale>0,'Invalid estimated similarity scale.');
translation=b-scale*a*Rrow;
transform=@(P) scale*double(P)*Rrow + translation;
newTie0=transform(oldTie);
newCam0=transform(oldCam);
errors=transform(oldMatch)-newMatch;
rms3D=sqrt(mean(sum(errors.^2,2)));
matchedNames=controlNames(matched);
checkTable=table(matchedNames,vecnorm(errors,2,2), errors(:,1),errors(:,2),errors(:,3),...
    'VariableNames',{'point_id','error3D_mm','dX_mm','dY_mm','dZ_mm'});
disp(checkTable);
fprintf('\nMatched points: %d; similarity scale = %.9f; fit RMS = %.3f mm\n', ...
    nnz(matched),scale,rms3D);
fprintf('WARNING: large residuals suggest wrong point IDs, poor target picking, or mismatched network geometry.\n');

% Object-to-camera rotation convention used here:
% R_old = Rz(kappa)*Ry(phi)*Rx(omega).
% For X_tls = s*Rcolumn*X_old + t, R_new = R_old*Rcolumn'.
% Verify THIS convention matches course full_matlab_bundle BEFORE use.
Rcolumn=Rrow';
newWpk=zeros(size(oldWpk));
for j=1:size(oldWpk,1)
    r=rotZ(oldWpk(j,3))*rotY(oldWpk(j,2))*rotX(oldWpk(j,1));
    rNew=r*Rcolumn';
    phi=atan2d(-rNew(3,1),hypot(rNew(1,1),rNew(2,1)));
    omega=atan2d(rNew(3,2),rNew(3,3));
    kappa=atan2d(rNew(2,1),rNew(1,1));
    newWpk(j,:)=[omega phi kappa];
end

if isfield(S,'activeImages')
    imageID=double(S.activeImages(:));
elseif isfield(S,'ext') && size(S.ext,1)==size(oldCam,1)
    imageID=double(S.ext(:,1));
else
    imageID=(1:size(oldCam,1))';
    warning('Image IDs not in MAT; assuming sequential IDs 1..N.');
end
assert(numel(imageID)==size(oldCam,1),'Camera count / image ID mismatch.');

outDir=uigetdir(fd,'Choose folder for Lab 5 initial-coordinate exports');
if isequal(outDir,0),return;end
prefix=fullfile(outDir,'lab5_tls_initial');
pointTable=table(names,newTie0(:,1),newTie0(:,2),newTie0(:,3), ...
    'VariableNames',{'point_id','X_mm','Y_mm','Z_mm'});
writetable(pointTable,[prefix '_tie_coordinates.csv']);
cameraTable=table(imageID,newCam0(:,1),newCam0(:,2),newCam0(:,3), ...
    newWpk(:,1),newWpk(:,2),newWpk(:,3), ...
    'VariableNames',{'image_id','Xc_mm','Yc_mm','Zc_mm','omega_deg','phi_deg','kappa_deg'});
writetable(cameraTable,[prefix '_camera_EOPs.csv']);
writetable(checkTable,[prefix '_similarity_residuals.csv']);

% Course format: image camera Xc Yc Zc omega phi kappa
ext=[imageID,99*ones(size(imageID)),newCam0,newWpk];
writeNumeric([prefix '.ext'],ext,'%d %d %.8f %.8f %.8f %.8f %.8f %.8f\n');
% Course bundle point numbers must come from Lab 4 results, NOT row indices.
assert(isfield(S,'numericID'), 'MAT lacks numericID used to write the .tie file.');
ids=double(S.numericID(:));
assert(numel(ids)==numel(names),'numericID/name mismatch.');
writeNumeric([prefix '.tie'],[ids,newTie0],'%d %.8f %.8f %.8f\n');
controlTable=table(controlNames(matched),newMatch(:,1),newMatch(:,2),newMatch(:,3), ...
    'VariableNames',{'point_id','X_mm','Y_mm','Z_mm'});
[~,controlIDrow]=ismember(controlNames(matched),names);
controlTable.BundlePointID=ids(controlIDrow);
writetable(controlTable,[prefix '_all_TLS_controls.csv']);

% Default THREE fixed controls selected to span 3D space; validate geometry.
chosen=["F1";"F7";"R12"];
assert(all(ismember(chosen,controlTable.point_id)), ...
    'Edit chosen control IDs in script to be among your CSV and Lab 4 network.');
[~,idx]=ismember(chosen,controlTable.point_id);
fixed=controlTable(idx,:);
writetable(fixed,[prefix '_three_FIXED_controls.csv']);
fixedXYZ=[fixed.X_mm,fixed.Y_mm,fixed.Z_mm];
triArea=norm(cross(fixedXYZ(2,:)-fixedXYZ(1,:), ...
    fixedXYZ(3,:)-fixedXYZ(1,:)))/2;
fprintf('Three-control triangle area = %.2f mm^2\n',triArea);
if triArea<1e4,warning('Control geometry is nearly collinear; select better points.');end

save([prefix '_transformation.mat'],'scale','Rrow','translation','newTie0', ...
    'newCam0','newWpk','names','ids','controlTable','checkTable','chosen','rms3D');
fprintf('\nSaved TLS initial EOPs, tie points, and candidate control points in %s\n',outDir);
fprintf('NEXT: verify omega/phi/kappa convention, generate course-compatible .con,\n');
fprintf('disable EOP priors, omit .dis, and run the LAB 5 bundle adjustment.\n');
end

function R=rotX(a)
c=cosd(a);s=sind(a);R=[1 0 0;0 c -s;0 s c];
end
function R=rotY(a)
c=cosd(a);s=sind(a);R=[c 0 s;0 1 0;-s 0 c];
end
function R=rotZ(a)
c=cosd(a);s=sind(a);R=[c -s 0;s c 0;0 0 1];
end
function writeNumeric(path,M,fmt)
f=fopen(path,'w');assert(f>=0,'Cannot write %s',path);
cleanup=onCleanup(@()fclose(f));
for k=1:size(M,1),fprintf(f,fmt,M(k,:));end
end
