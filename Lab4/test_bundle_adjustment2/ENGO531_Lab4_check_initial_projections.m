function ENGO531_Lab4_check_initial_projections()
% ENGO 531 Lab 4: compare initial right-wall XYZ with image measurements.
% Choose the *_results.mat files made by ENGO531_Lab4_resection_warmstart_FIXED.
% Independent projection diagnostic: NO bundle adjustment, NO file editing.
clc;
[f6,p6]=uigetfile('*_results.mat','Select SIX-image result MAT (warm_images_1_2_3_4_5_6_results.mat)');
if isequal(f6,0),return;end
[f7,p7]=uigetfile('*_results.mat','Select SEVEN-image result MAT (warm_images_1_2_3_4_5_6_7_results.mat)',p6);
if isequal(f7,0),return;end
A=load(fullfile(p6,f6)); B=load(fullfile(p7,f7));
assert(isfield(A,'eop') && isfield(A,'tie') && isfield(A,'names'), ...
    'Six-image MAT must contain eop, tie, names');
assert(isfield(B,'XYZ0') && isfield(B,'names') && isfield(B,'pointName') && ...
    isfield(B,'imageID') && isfield(B,'xy') && isfield(B,'ext'), ...
    'Seven-image MAT missing initialization and observation data');
assert(size(A.eop.XYZc,1)==6 && size(B.ext,1)==7, ...
    'Please select SIX-image and SEVEN-image MAT files in that order.');
IOP = B.IOP(:)';
% The custom resection used the following rotation construction, but the
% course P-code's collinearity convention is not yet independently verified.
% We test two signs and two transpose directions against adjusted F points.
rotForms = {'RzT_Ry_RxT','Rz_Ry_Rx','Rz_RyT_RxT','RzT_RyT_Rx'};
ySigns=[-1,1];
image6pose=[A.eop.XYZc(6,:),A.eop.wpk(6,:)];
is6=B.imageID==6;
baseIsF=startsWith(B.pointName,'F') & is6;
refNames=B.pointName(baseIsF); refObs=B.xy(baseIsF,:);
refX=zeros(numel(refNames),3); hasRef=false(numel(refNames),1);
for j=1:numel(refNames)
    idx=find(A.names==refNames(j),1);
    if ~isempty(idx)
        refX(j,:)=A.tie.XYZ(idx,:);hasRef(j)=true;
    end
end
refNames=refNames(hasRef);refObs=refObs(hasRef,:);refX=refX(hasRef,:);
assert(numel(refNames)>=4,'Insufficient adjusted front-face points in Image 6 to verify projection.');
fprintf('\n=== VERIFY PROJECTION AGAINST ADJUSTED F TARGETS IN IMAGE 6 ===\n');
results=table('Size',[numel(rotForms)*2 4],'VariableTypes', ...
    {'string','double','double','double'}, ...
    'VariableNames',{'Rotation','YSign','RMS_px','PositiveDepthCount'});
q=0;
for i=1:numel(rotForms)
    for sy=ySigns
        q=q+1;
        [pred,depth]=project(refX,image6pose,IOP,rotForms{i},sy);
        err=sqrt(sum((pred-refObs).^2,2));
        results(q,:)={string(rotForms{i}),sy,sqrt(mean(err.^2)),sum(depth>0)};
    end
end
results=sortrows(results,'RMS_px');disp(results);
form=char(results.Rotation(1));sy=results.YSign(1);
fprintf('Best of the candidate conventions: %s; Y sign %d.\n',form,sy);
fprintf('Calibration RMS on adjusted front points: %.3f px\n',results.RMS_px(1));
if results.RMS_px(1)>15
    warning(['Projection convention is NOT VERIFIED: even adjusted front points miss by %.1f px. ' ...
    'DO NOT interpret right-face projection errors as bad XYZ until the course collinearity convention is established.'],results.RMS_px(1));
end

% Report every retained right-face point independently for Images 6 and 7.
% Image 6 uses adjusted camera EOPs. Image 7 uses field-note EOPs and the
% approximate/resected EOP that the previous script wrote into B.ext.
right=startsWith(B.names,'R'); rn=B.names(right); rXYZ=B.XYZ0(right,:);
idx6=find(B.ext(:,1)==6,1);idx7=find(B.ext(:,1)==7,1);
assert(~isempty(idx6)&&~isempty(idx7),'Seven-image MAT lacks Image 6 or Image 7');
pose7notes=[4368,-97,-495,90,60,0];
pose7resected=B.ext(idx7,3:8);
T6=makeReport(rn,rXYZ,B,6,image6pose,IOP,form,sy);
T7initial=makeReport(rn,rXYZ,B,7,pose7notes,IOP,form,sy);
T7resected=makeReport(rn,rXYZ,B,7,pose7resected,IOP,form,sy);

fprintf('\n=== IMAGE 6: converged EOP + APPROXIMATE R XYZ ===\n');
disp(sortrows(T6,'Error_px','descend'));
fprintf('\n=== IMAGE 7: FIELD-NOTE EOP + APPROXIMATE R XYZ ===\n');
disp(sortrows(T7initial,'Error_px','descend'));
fprintf('\n=== IMAGE 7: PREVIOUS RESECTION EOP + APPROXIMATE R XYZ ===\n');
disp(sortrows(T7resected,'Error_px','descend'));

outDir=p7;
writetable(results,fullfile(outDir,'projection_convention_check.csv'));
writetable(T6,fullfile(outDir,'projection_R_image6.csv'));
writetable(T7initial,fullfile(outDir,'projection_R_image7_fieldnotes.csv'));
writetable(T7resected,fullfile(outDir,'projection_R_image7_resected.csv'));
fprintf('\nSaved four diagnostic CSV files to: %s\n',outDir);
fprintf(['If reference F point error is low but Image 6 R errors are high, ' ...
    'recheck R-wall XYZ initialization.\n']);
fprintf(['If reference F and Image 6 R errors are low, focus on Image 7 EOP ' ...
    'and Image 7 measurements.\n']);
fprintf('A poor projection check does not prove any single observed point is wrong.\n');
end

function T=makeReport(names,X,B,imageNum,pose,IOP,form,ysign)
[pred,depth]=project(X,pose,IOP,form,ysign);
obs=nan(numel(names),2);
for j=1:numel(names)
    k=find(B.pointName==names(j) & B.imageID==imageNum,1);
    if ~isempty(k), obs(j,:)=B.xy(k,:);end
end
err=hypot(pred(:,1)-obs(:,1),pred(:,2)-obs(:,2));
T=table(names,obs(:,1),obs(:,2),pred(:,1),pred(:,2),err,depth,depth>0, ...
 'VariableNames',{'Target','ObservedX_px','ObservedY_px','PredictedX_px', ...
 'PredictedY_px','Error_px','CameraDepth','PositiveDepth'});
T=T(all(isfinite(obs),2),:);
end

function [pixel,depth]=project(X,pose,IOP,which,ysign)
w=deg2rad(pose(4));p=deg2rad(pose(5));k=deg2rad(pose(6));
Rx=[1 0 0;0 cos(w) -sin(w);0 sin(w) cos(w)];
Ry=[cos(p) 0 sin(p);0 1 0;-sin(p) 0 cos(p)];
Rz=[cos(k) -sin(k) 0;sin(k) cos(k) 0;0 0 1];
switch which
    case 'RzT_Ry_RxT', M=Rz'*Ry*Rx';
    case 'Rz_Ry_Rx',   M=Rz*Ry*Rx;
    case 'Rz_RyT_RxT', M=Rz*Ry'*Rx';
    case 'RzT_RyT_Rx', M=Rz'*Ry'*Rx;
    otherwise, error('Unknown rotation option');
end
Q=M*(X-pose(1:3))'; depth=Q(3,:)';
d=depth;d(abs(d)<1e-9)=NaN;
pixel=[IOP(1)-IOP(3)*Q(1,:)'./d, IOP(2)+ysign*IOP(3)*Q(2,:)'./d];
end
