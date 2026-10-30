function ENGO531_Lab4_final_report_NO_F1_F6()
% ENGO 531 Lab 4: post-adjustment analysis and independent distance checks.
% Run AFTER bundle adjustment/data snooping. Does not change observations.
% Supports *_results.mat from the cumulative snooping script and joint run.
clc;
[file,folder]=uigetfile('*.mat','Select FINAL bundle adjustment *_results.mat');
if isequal(file,0),return;end
S=load(fullfile(folder,file));
req={'eop','tie','Cx','vr_img'};
assert(all(isfield(S,req)),'MAT file needs eop, tie, Cx and vr_img');
[~,base,~]=fileparts(file);
base=regexprep(base,'_results$','');
output=fullfile(folder,[base '_REPORT']);
if ~isfolder(output),mkdir(output);end
XYZ=double(S.tie.XYZ); CAM=double(S.eop.XYZc);
ang=double(S.eop.wpk); nP=size(XYZ,1); nC=size(CAM,1);
assert(size(XYZ,2)==3 && size(CAM,2)==3,'XYZ arrays must have three columns');
% Point identifiers: use MAT names when available, otherwise map CSV.
if isfield(S,'names') && numel(S.names)==nP
    labels=string(S.names(:));
else
    labels=strings(nP,1);
    mapFile=fullfile(folder,[base '_target_ID_map.csv']);
    if ~isfile(mapFile),mapFile=fullfile(folder,[base '_point_map.csv']);end
    if isfile(mapFile)
       M=readtable(mapFile,'TextType','string');
       if ismember('OriginalPointID',M.Properties.VariableNames)
          mapNames=string(M.OriginalPointID); mapIDs=string(M.BundlePointID);
       elseif all(ismember({'Target','ID'},M.Properties.VariableNames))
          mapNames=string(M.Target);mapIDs=string(M.ID);
       else
          mapNames=strings(0,1);mapIDs=strings(0,1);
       end
       if isfield(S.tie,'pt')
           theseIDs=string(S.tie.pt(:));
       elseif isfield(S,'numericID')
           theseIDs=string(S.numericID(:));
       else
           theseIDs=string((1001:1000+nP)');
       end
       [found,loc]=ismember(theseIDs,mapIDs);
       labels(found)=mapNames(loc(found));
    end
    missing=labels=="";
    if isfield(S.tie,'pt')
        tmp=string(S.tie.pt(:));labels(missing)=tmp(missing);
    else
        labels(missing)="P"+string(find(missing));
    end
end
assert(numel(unique(labels))==nP,'Duplicate target names: check mapping and tie ordering');
imageIDs=(1:nC)';
if isfield(S.eop,'img')
    parsed=str2double(string(S.eop.img(:)));
    if all(isfinite(parsed)) && numel(parsed)==nC,imageIDs=parsed;end
end
% Observations from the FINAL pass, not pre-snooping PHO.
if isfield(S,'photoTable'), PHO=double(S.photoTable);
else
    src=fullfile(folder,[base '.pho']);
    if ~isfile(src)
        [pf,pd]=uigetfile({'*.pho;*.PHO','Final numeric bundle PHO'},'Select FINAL .pho file',folder);
        assert(~isequal(pf,0),'Need final PHO file to label image residuals');src=fullfile(pd,pf);
    end
    PHO=readmatrix(src,'FileType','text');
end
assert(size(PHO,2)>=4,'PHO must have [pointID,imageID,x,y]');
v=double(S.vr_img(:,1));
assert(numel(v)==2*size(PHO,1),'Residual/PHO row count differs: select matching FINAL PHO');
V=reshape(v,2,[])';imageRMS=sqrt(mean(v.^2));
% Observation sigmas: reconstruct from final CFG if present.
if isfield(S,'cfg'),cfg=double(S.cfg(:));
else
    cf=fullfile(folder,[base '.cfg']);
    if isfile(cf),cfg=readmatrix(cf,'FileType','text');cfg=cfg(:);
    else,cfg=[1;0;0.001;0.0001;0.01;3.2;5;1;0.01];
        warning('CFG missing: using fallback SDs [3.2 px,5 mm,1 mm,0.01 deg]');
    end
end
assert(numel(cfg)>=9,'CFG is incomplete');
getv=@(field) getResidual(S,field);
ve=getv('vr_eop');vd=getv('vr_dis');vc=getv('vr_con');
% Explicit weighted residual sums; control residuals require separate SD.
vPvImage=sum((v./cfg(6)).^2);
sigE=repmat([cfg(8);cfg(8);cfg(8);cfg(9);cfg(9);cfg(9)],ceil(numel(ve)/6),1);
vPvEOP=sum((ve./sigE(1:numel(ve))).^2);
vPvDist=sum((vd./cfg(7)).^2);
if ~isempty(vc)
    warning('Control residuals exist but their SD is unspecified; variance factor omitted.');
    validFactor=false;vPv=NaN;
else
    validFactor=true;vPv=vPvImage+vPvEOP+vPvDist;
end
nObs=numel(v)+numel(ve)+numel(vd)+numel(vc);
nUnknown=size(S.Cx,1);DOF=nObs-nUnknown;
if DOF<=0,validFactor=false;warning('Nonpositive nominal DOF');end
if validFactor,sigma0sq=vPv/DOF;sigma0=sqrt(sigma0sq);
else,sigma0sq=NaN;sigma0=NaN;end
% Cumulative factor is against the original SDs; local factor uses final CFG.
if isfield(S,'finalCumulativeVarianceScale')
  cumulative=double(S.finalCumulativeVarianceScale);
elseif isfield(S,'updatedCumulativeVarianceScale')
  cumulative=double(S.updatedCumulativeVarianceScale);
elseif isfield(S,'scaleHistory') && ~isempty(S.scaleHistory)
  cumulative=double(S.scaleHistory(end,6));
else
  cumulative=NaN;
end
fprintf('\nFINAL BUNDLE REPORT\nImages %d; targets %d; image points %d\n',nC,nP,size(PHO,1));
fprintf('Image RMS per scalar coordinate: %.4f px\n',imageRMS);
fprintf('n=%d, u=%d, DOF=%d, weighted vPv=%.6f\n',nObs,nUnknown,DOF,vPv);
fprintf('Local a posteriori variance factor: %.6f; sigma0: %.6f\n',sigma0sq,sigma0);
if isfinite(cumulative),fprintf('Recorded cumulative variance factor: %.6f\n',cumulative);end
summary=table(nC,nP,size(PHO,1),nObs,nUnknown,DOF,imageRMS,vPv, ...
    vPvImage,vPvEOP,vPvDist,sigma0sq,sigma0,cumulative, ...
    'VariableNames',{'Images','Targets','PhotoPoints','ScalarObs','Unknowns','DOF', ...
    'ImageRMS_px','vPv','vPvImage','vPvEOP','vPvDistance','LocalVarFactor','LocalSigma0','CumulativeVarFactor'});
disp(summary);writetable(summary,fullfile(output,'Adjustment_summary.csv'));
% Standardized image residuals, with CURRENT sigma_xy from final CFG.
if size(S.vr_img,2)>=2,red=double(S.vr_img(:,2));
else,red=nan(size(v));end
w=nan(size(v));ok=red>0 & isfinite(red);w(ok)=v(ok)./(cfg(6)*sqrt(red(ok)));
W=reshape(w,2,[])';R=reshape(red,2,[])';
% Make names from point ID, but do not assume IDs match row order.
if isfield(S.tie,'pt'), ptids=string(S.tie.pt(:));
elseif isfield(S,'numericID'),ptids=string(S.numericID(:));
else,ptids=string((1001:1000+nP)');end
[found,loc]=ismember(string(PHO(:,1)),ptids);
obsLabels=strings(size(PHO,1),1);obsLabels(found)=labels(loc(found));obsLabels(~found)=string(PHO(~found,1));
maxW=max(abs(W),[],2);flag=maxW>3;
residualTable=table(obsLabels,PHO(:,1),PHO(:,2),PHO(:,3),PHO(:,4), ...
    V(:,1),V(:,2),R(:,1),R(:,2),W(:,1),W(:,2),maxW,flag, ...
    'VariableNames',{'Target','ID','Image','x_px','y_px','vx_px','vy_px', ...
    'rX','rY','wX','wY','MaxAbsW','FlagAbsWgt3'});
writetable(residualTable,fullfile(output,'Image_residuals.csv'));
disp('Largest standardized residuals:');
sorted=sortrows(residualTable,'MaxAbsW','descend');disp(sorted(1:min(15,height(sorted)),:));
ims=unique(PHO(:,2));ir=zeros(numel(ims),2);
for k=1:numel(ims)
    mask=PHO(:,2)==ims(k);vv=V(mask,:);
    ir(k,:)=[nnz(mask),sqrt(mean(vv(:).^2))];
end
perImage=table(ims,ir(:,1),ir(:,2),'VariableNames',{'Image','PhotoPoints','RMS_px'});
writetable(perImage,fullfile(output,'Image_RMS.csv'));disp(perImage);
% Coordinates and precision (assuming EOP6, then tie3 indexing).
C=double(S.Cx);hasCov=isequal(size(C),[6*nC+3*nP,6*nC+3*nP]);
if ~hasCov,warning('Unexpected Cx dimensions: SDs/ellipses disabled');end
sdP=nan(nP,3);sdC=nan(nC,6);
if hasCov
 sdAll=sqrt(max(diag(C),0));sdC=reshape(sdAll(1:6*nC),6,[])';
 sdP=reshape(sdAll(6*nC+1:end),3,[])';
end
T=table(labels,XYZ(:,1),XYZ(:,2),XYZ(:,3),sdP(:,1),sdP(:,2),sdP(:,3), ...
  'VariableNames',{'Target','X_mm','Y_mm','Z_mm','sX_mm','sY_mm','sZ_mm'});
writetable(T,fullfile(output,'Adjusted_points.csv'));
eopT=array2table([imageIDs CAM ang sdC], 'VariableNames', ...
 {'Image','Xc_mm','Yc_mm','Zc_mm','Omega_deg','Phi_deg','Kappa_deg', ...
 'sXc_mm','sYc_mm','sZc_mm','sOmega_deg','sPhi_deg','sKappa_deg'});
writetable(eopT,fullfile(output,'Adjusted_EOPs.csv'));
% PRESET FIELD DISTANCES (metres converted to millimetres).
% The list is recreated every run; no prompts or stale input CSVs.
% F1-F6 intentionally excluded from comparison (used as a scale observation).
from = ["L9"; "L2"; "R10"; "R4"; "F3"];
to   = ["L10"; "L4"; "R9"; "R2"; "F4"];
tape_m = [0.816; 1.798; 0.816; 1.798; 0.967];
tape = 1000*tape_m;
D=table(from,to,tape_m,tape,'VariableNames',{'From','To','Tape_m','Tape_mm'});
checkFile=fullfile(output,'Check_distances_input.csv');
writetable(D,checkFile);
fprintf('\nComparing %d supplied measured distances (input metres, calculations mm).\n',height(D));
[ok1,a]=ismember(string(D.From),labels);
[ok2,b]=ismember(string(D.To),labels);
if any(~ok1|~ok2)
    disp('Pairs with point labels absent from the FINAL adjustment:');
    disp(D(~ok1|~ok2,:));
    warning('Skipping check distances with absent endpoints.');
    D=D(ok1&ok2,:); a=a(ok1&ok2); b=b(ok1&ok2);
end
if height(D)>0
    % Identify the two measured distances already present in the adjustment.
    % Always treat F1-F6 and R4-R1 as used, even if the .dis file is missing.
    usedInBundle = (D.From=="F1" & D.To=="F6") | ...
                   (D.From=="F6" & D.To=="F1") | ...
                   (D.From=="R4" & D.To=="R1") | ...
                   (D.From=="R1" & D.To=="R4");
    usedRows=[];
    if isfield(S,'distanceRows'),usedRows=double(S.distanceRows);
    else
        disPath=fullfile(folder,[base '.dis']);
        if isfile(disPath),usedRows=readmatrix(disPath,'FileType','text');end
    end
    if ~isempty(usedRows) && size(usedRows,2)>=2
        usedPairs=sort(string(usedRows(:,1:2)),2);
        for jj=1:height(D)
            pair=sort([ptids(a(jj)) ptids(b(jj))]);
            if any(all(usedPairs==pair,2))
                usedInBundle(jj)=true;
            end
        end
    end
    measured=double(D.Tape_mm);
    adjusted=nan(height(D),1);delta=adjusted;relppm=adjusted;sdDist=adjusted;
    for k=1:height(D)
        diffXYZ=XYZ(b(k),:)-XYZ(a(k),:);L=norm(diffXYZ);
        adjusted(k)=L;delta(k)=L-measured(k);relppm(k)=1e6*delta(k)/measured(k);
        if hasCov && L>0
            ia=6*nC+(a(k)-1)*3+(1:3);ib=6*nC+(b(k)-1)*3+(1:3);
            unit=diffXYZ(:)/L;
            Cdiff=C(ia,ia)+C(ib,ib)-C(ia,ib)-C(ib,ia);
            sdDist(k)=sqrt(max(0,unit'*Cdiff*unit));
        end
    end
    checkType=repmat("Independent",height(D),1);
    checkType(usedInBundle)="Used in bundle";
    checks=table(string(D.From),string(D.To),double(D.Tape_m),measured,adjusted, ...
        delta,abs(delta),relppm,sdDist,checkType, ...
        'VariableNames',{'From','To','Tape_m','Tape_mm','Adjusted_mm', ...
        'AdjustedMinusTape_mm','AbsoluteError_mm','Error_ppm', ...
        'PropagatedSD_mm','CheckType'});
    disp('DISTANCE COMPARISON (mm; used-in-bundle distances are NOT independent):');
    disp(checks);
    writetable(checks,fullfile(output,'Check_distance_comparison.csv'));
    independent = ~usedInBundle;
    if any(independent)
        fprintf('Independent check-distance RMSE = %.3f mm (n=%d)\n', ...
            sqrt(mean(delta(independent).^2)),nnz(independent));
        fprintf('Independent check-distance mean error = %.3f mm\n', ...
            mean(delta(independent)));
    else
        warning('No independent check distances available to calculate RMSE.');
    end
    f=figure('Name','Distance differences');
    namesForPlot=string(D.From)+"-"+string(D.To);
    bar(categorical(namesForPlot),delta);
    ylabel('Adjusted minus tape (mm)');grid on;
    title('All field distances (see CSV for bundle vs independent classification)');
    saveas(f,fullfile(output,'Check_distances.png'));
else
    warning('No supplied distance pairs had both endpoints in the final network.');
end
% NETWORK PLOTS
colors=[0.1 .45 .85;.12 .65 .38;.9 .45 .16];prefixes=["L","F","R"];
f=figure('Name','Final adjusted network');hold on;grid on;axis equal;view(3);
for k=1:3
 mask=startsWith(labels,prefixes(k));
 if any(mask)
 scatter3(XYZ(mask,1),XYZ(mask,2),XYZ(mask,3),35,colors(k,:),'filled');
 end
end
text(XYZ(:,1),XYZ(:,2),XYZ(:,3),cellstr(labels),'FontSize',8);
plot3(CAM(:,1),CAM(:,2),CAM(:,3),'k--^','LineWidth',1);
for k=1:nC,text(CAM(k,1),CAM(k,2),CAM(k,3),sprintf(' Camera %d',imageIDs(k)));end
xlabel('X (mm)');ylabel('Y (mm)');zlabel('Z (mm)');title('Final bundle: object points and cameras');
saveas(f,fullfile(output,'Network_3D.png'));
f=figure('Name','Image residuals');tiledlayout(3,ceil(numel(ims)/3));
for k=1:numel(ims)
 nexttile;mask=PHO(:,2)==ims(k);hold on;grid on;
 scatter(PHO(mask,3),PHO(mask,4),10,'k','filled');
 quiver(PHO(mask,3),PHO(mask,4),30*V(mask,1),30*V(mask,2),0,'r');
 title(sprintf('Image %d',ims(k)));set(gca,'YDir','reverse');axis equal;
 xlim([0 5712]);ylim([0 4284]);
end
sgtitle('Final image residuals (30x vectors)');saveas(f,fullfile(output,'Residual_vectors.png'));
f=figure('Name','Standardized residuals');plot(w,'.');hold on;grid on;
yline(3,'r--');yline(-3,'r--');xlabel('Scalar observation index');ylabel('w');
title('Final standardized image residuals');saveas(f,fullfile(output,'Standardized_residuals.png'));
% Approximate 95% confidence regions assuming Cx is an absolute covariance matrix.
% If Cx is a cofactor matrix, uncertainties require sigma0 scaling instead.
if hasCov
 f=figure('Name','95% point XY ellipses');hold on;grid on;axis equal;
 scatter(XYZ(:,1),XYZ(:,2),15,'k','filled');
 for k=1:nP
    ix=6*nC+(k-1)*3+(1:2);
    plotEllipse(C(ix,ix),XYZ(k,1:2),sqrt(5.991464547),3); %scale ellipse
 end
 text(XYZ(:,1),XYZ(:,2),cellstr(labels),'FontSize',8);
 xlabel('X (mm)');ylabel('Y (mm)');title('95% point XY confidence ellipses (10x display)');
 saveas(f,fullfile(output,'Point_XY_ellipses.png'));
 f=figure('Name','95% camera XYZ ellipsoids');hold on;grid on;axis equal;view(3);
 for k=1:nC
    ix=(k-1)*6+(1:3);
    plotEllipsoid(C(ix,ix),CAM(k,:),sqrt(7.814727903),5);
    plot3(CAM(k,1),CAM(k,2),CAM(k,3),'k.');
 end
 xlabel('X (mm)');ylabel('Y (mm)');zlabel('Z (mm)');
 title('95% camera position ellipsoids (5x display)');
 saveas(f,fullfile(output,'Camera_XYZ_ellipsoids.png'));
end
fprintf('\nSaved reports and plots to: %s\n',output);
fprintf('NOTE: convergence must have been explicitly confirmed in the bundle output.\n');
fprintf('NOTE: check Cx covariance convention before interpreting precision/ellipses.\n');
end
function v=getResidual(S,field)
 if ~isfield(S,field)||isempty(S.(field)),v=[];return;end
 q=S.(field);if ~isnumeric(q) || (isscalar(q) && q==0),v=[];return;end
 v=double(q(:,1));
end
function plotEllipse(C,center,confidence,displayScale)
 C=(C+C')/2;[V,D]=eig(C);d=diag(D);
 if any(~isfinite(d))||any(d<-1e-8),return;end
 t=linspace(0,2*pi,80);
 P=V*diag(sqrt(max(d,0))*confidence*displayScale)*[cos(t);sin(t)];
 plot(center(1)+P(1,:),center(2)+P(2,:),'Color',[.7 .3 .2]);
end
function plotEllipsoid(C,center,confidence,displayScale)
 C=(C+C')/2;[V,D]=eig(C);d=diag(D);
 if any(~isfinite(d))||any(d<-1e-8),return;end
 [x,y,z]=sphere(16);A=V*diag(sqrt(max(d,0))*confidence*displayScale);
 pts=A*[x(:)';y(:)';z(:)'];
 surf(reshape(pts(1,:),size(x))+center(1), ...
      reshape(pts(2,:),size(y))+center(2), ...
      reshape(pts(3,:),size(z))+center(3), ...
      'FaceAlpha',0.15,'EdgeColor','none');
end
