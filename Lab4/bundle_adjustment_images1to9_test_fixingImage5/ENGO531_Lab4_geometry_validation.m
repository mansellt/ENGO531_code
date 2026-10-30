function ENGO531_Lab4_geometry_validation()
% ENGO 531 Lab 4: validate collinearity before diagnosing the right wall.
% Reads existing six-image (CONVERGED) and seven-image result MAT files.
% Does not change input measurements or run a new adjustment.
clc;
[f6,d6]=uigetfile('*_results.mat','Select CONVERGED six-image _results.mat');
if isequal(f6,0),return;end
[f7,d7]=uigetfile('*_results.mat','Select seven-image _results.mat (initial XYZ)',d6);
if isequal(f7,0),return;end
A=load(fullfile(d6,f6)); B=load(fullfile(d7,f7));
assert(all(isfield(A,{'eop','tie','names'})),'Six-image file missing eop/tie/names');
assert(all(isfield(B,{'names','XYZ0','ext','IOP'})),'Seven-image file missing names/XYZ0/ext/IOP');
assert(size(A.eop.XYZc,1)==6 && size(B.ext,1)==7,'Wrong six/seven image files');
% Use numeric photo/target IDs to avoid mismatching adjusted point ordering.
% Prefer saving names/numericID and photoTable in the source MAT as done by course scripts.
assert(isfield(A,'numericID') && isfield(A,'photoTable'), ...
 'Six-image MAT must contain numericID and photoTable. Regenerate from the supplied script.');
assert(isfield(B,'numericID') && isfield(B,'photoTable'), ...
 'Seven-image MAT must contain numericID and photoTable.');
[Xref,ids]=adjustedXYZ(A);
assert(size(Xref,1)==numel(ids),'Adjusted XYZ/ID size disagreement');
T6=A.photoTable; isRef=ismember(T6(:,1),ids) & T6(:,2)>=2 & T6(:,2)<=6;
T6=T6(isRef,:); [has,j]=ismember(T6(:,1),ids);assert(all(has));
XYZ=Xref(j,:); obs=T6(:,3:4); im=T6(:,2);
valid=all(isfinite(XYZ),2) & all(isfinite(obs),2);
XYZ=XYZ(valid,:);obs=obs(valid,:);im=im(valid);
fprintf('\nVALIDATION DATA: %d reference photo observations from images 2-6.\n',numel(im));
% IMPORTANT: EOP/order could be unusual. Search camera rotation ordering,
% transpose, Euler sign combinations, and camera optical axis choice.
% A projection convention must reproduce reference F/L points first.
sgn=[-1 1]; opts=[]; n=0;
for order=1:6
 for trans=0:1
  for sw=sgn
   for sp=sgn
    for sk=sgn
     for depthAxis=1:3
      for imageX=1:3
       if imageX==depthAxis,continue;end
       imageY=6-depthAxis-imageX;
       for sx=sgn
        for sy=sgn
         for sf=sgn
          n=n+1;
          cfg=[order trans sw sp sk depthAxis imageX imageY sx sy sf];
          [nGood,e]=evalErrors(XYZ,obs,im,A.eop,B.IOP,cfg);
          opts(n,:)=[cfg nGood e]; %#ok<AGROW>
         end
        end
       end
      end
     end
    end
   end
  end
 end
end
% RMS in last position; count of finite pixels second to last
opts=sortrows(opts,size(opts,2));
Top=array2table(opts(1:10,:), 'VariableNames', ...
 {'EulerOrder','Transpose','OmegaSign','PhiSign','KappaSign','DepthAxis', ...
 'ImageXAxis','ImageYAxis','XSign','YSign','FocalSign','ValidCount','RMS_px'});
fprintf('\nBEST CAMERA CONVENTIONS (against known adjusted tie points):\n');disp(Top);
writetable(Top,fullfile(d7,'projection_validation_candidates.csv'));
best=opts(1,1:11);bestRMS=opts(1,end);nValid=opts(1,end-1);
% Require almost all ref projections and roughly correct adjusted fit.
if nValid < .95*numel(im) || ~isfinite(bestRMS) || bestRMS>15
 fprintf(2,['\nSTOP: not validated (best RMS %.2f px; valid %d/%d).\n' ...
  'Likely adjusted tie-point ordering, EOP/angle convention, or image coordinate convention.\n' ...
  'Do not interpret right-wall reprojection errors from an unverified model.\n'], ...
  bestRMS,nValid,numel(im));
 return;
end
fprintf('\nProjection model VERIFIED on reference data: RMS %.3f px.\n',bestRMS);
% Compare front-right corner and wall geometry before target error checks.
frontNames=["F1" "F6" "F7" "F12"];
fprintf('\nADJUSTED FRONT POINTS (for coordinate-frame comparison):\n');
for k=1:numel(frontNames)
 j=find(A.names==frontNames(k),1);
 if ~isempty(j)
  fprintf('%-5s  %10.2f %10.2f %10.2f mm\n',frontNames(k),Xref(j,:));
 end
end
rMask=startsWith(B.names,'R');rNames=B.names(rMask);rXYZ=B.XYZ0(rMask,:);
fprintf('\nRIGHT-WALL CORNERS (approximate, not control):\n');
for nm=["R4" "R1" "R5" "R8"]
 j=find(rNames==nm,1);
 if ~isempty(j),fprintf('%-5s  %10.2f %10.2f %10.2f mm\n',nm,rXYZ(j,:));end
end
% Report measured-vs-predicted right points for 6/7; image 6 uses converged EOP.
for imageNum=[6 7]
 if imageNum==6
  pose=[A.eop.XYZc(6,:),A.eop.wpk(6,:)];
 else
  r=find(B.ext(:,1)==7,1);pose=B.ext(r,3:8);
 end
 [pred,depth]=projectPoints(rXYZ,pose,B.IOP,best);
 id= B.numericID(rMask); [m,j]=ismember(B.photoTable(:,1),id);
 rows=(B.photoTable(:,2)==imageNum)&m; P=B.photoTable(rows,:);
 [~,which]=ismember(P(:,1),id);
 errorXY=pred(which,:)-P(:,3:4);
 radial=hypot(errorXY(:,1),errorXY(:,2));
 report=table(rNames(which),P(:,3),P(:,4),pred(which,1),pred(which,2), ...
  radial,depth(which),'VariableNames',{'Target','ObservedX','ObservedY', ...
  'PredictedX','PredictedY','Error_px','DepthCoordinate'});
 report=sortrows(report,'Error_px','descend');
 fprintf('\nIMAGE %d: right wall initial-geometry reprojection test:\n',imageNum);disp(report);
 writetable(report,fullfile(d7,sprintf('verified_projection_R_image%d.csv',imageNum)));
end
fprintf('\nValidation passed; compare Image 6 versus Image 7 error patterns before altering any observations.\n');
end

function [X,ids]=adjustedXYZ(A)
X=A.tie.XYZ; ids=A.numericID(:);
% If returned tie output includes IDs, prefer actual returned order.
if isfield(A.tie,'id'),ids=A.tie.id(:);end
if isfield(A.tie,'ID'),ids=A.tie.ID(:);end
if numel(ids)~=size(X,1)
 error('Could not establish matching tie IDs and adjusted XYZ order.');
end
end

function [nvalid,rms]=evalErrors(X,obs,im,E,IOP,cfg)
res=zeros(size(obs));ok=true(size(im));
for k=unique(im(:))'
 p=[E.XYZc(k,:),E.wpk(k,:)];i=im==k;
 [pred,~]=projectPoints(X(i,:),p,IOP,cfg);
 res(i,:)=pred-obs(i,:);ok(i)=all(isfinite(pred),2);
end
nvalid=sum(ok);if nvalid==0,rms=Inf;else
 e=hypot(res(ok,1),res(ok,2));rms=sqrt(mean(e.^2));end
if nvalid<.95*numel(im),rms=Inf;end
end

function [pix,depth]=projectPoints(X,pose,IOP,cfg)
angles=deg2rad(pose(4:6)).*cfg(3:5);
a=angles(1);b=angles(2);c=angles(3);
Rx=[1 0 0;0 cos(a) -sin(a);0 sin(a) cos(a)];
Ry=[cos(b) 0 sin(b);0 1 0;-sin(b) 0 cos(b)];
Rz=[cos(c) -sin(c) 0;sin(c) cos(c) 0;0 0 1];
R={Rx,Ry,Rz};v=perms(1:3);v=v(cfg(1),:);
M=R{v(1)}*R{v(2)}*R{v(3)};
if cfg(2)==1,M=M';end
Q=(M*(X-pose(1:3))')';
depth=Q(:,cfg(6));depth(abs(depth)<1e-8)=NaN;
pix=[IOP(1)+cfg(11)*cfg(9)*IOP(3)*Q(:,cfg(7))./depth, ...
     IOP(2)+cfg(11)*cfg(10)*IOP(3)*Q(:,cfg(8))./depth];
end
