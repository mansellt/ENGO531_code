function ENGO531_Lab5_compare_TLS_checkpoints()
% ENGO 531 Lab 5: independent TLS check-point comparison of TWO saved
% MODE-0 bundle adjustments (e.g., three controls versus six controls).
%
% Inputs:
%   1) Six-control lab5_TLS_results.mat
%   2) Three-control lab5_TLS_results.mat
%   3) target_coordinates(3).csv [point_id,X,Y,Z,cloud_index]
%
% IMPORTANT:
% - TLS CSV X,Y,Z are METRES; bundle XYZ coordinates are MILLIMETRES.
% - A TLS point is an INDEPENDENT check ONLY if it was NOT used as a
%   control in that solution. The code excludes each solution's controls.
% - Fair head-to-head results use ONLY points that are checks in BOTH runs.
% - The initial .tie coordinates are NOT used: adjusted tie.XYZ is used.
% - Check points are independent of the control constraints, but are still
%   image observations participating in the bundle adjustment.
% - This code does not rerun or change either adjustment.

clc;
fprintf('\nENGO 531 Lab 5: independent TLS check-point comparison\n');
[file3,dir3] = uigetfile('*.mat','Select THREE-control Lab 5 results MAT');
if isequal(file3,0), return; end
[file6,dir6] = uigetfile('*.mat','Select SIX-control Lab 5 results MAT');
if isequal(file6,0), return; end
[csvFile,csvDir] = uigetfile('*.csv','Select TLS target_coordinates CSV');
if isequal(csvFile,0), return; end
outDir = uigetdir(pwd,'Choose folder for comparison CSVs and figures');
if isequal(outDir,0), return; end

TLS = readtable(fullfile(csvDir,csvFile),'TextType','string');
needed = {'point_id','X','Y','Z'};
assert(all(ismember(needed,TLS.Properties.VariableNames)), ...
    'TLS CSV needs point_id,X,Y,Z (XYZ in metres).');
TLS.Target = upper(strtrim(string(TLS.point_id)));
assert(numel(unique(TLS.Target))==height(TLS), 'Duplicate TLS target names.');
tlsXYZ = 1000 * double([TLS.X,TLS.Y,TLS.Z]);  % m -> mm
assert(all(isfinite(tlsXYZ(:))), 'TLS CSV has invalid XYZ values.');

A = extractChecks(fullfile(dir3,file3),TLS.Target,tlsXYZ,'3 controls');
B = extractChecks(fullfile(dir6,file6),TLS.Target,tlsXYZ,'6 controls');

writetable(A.Checks,fullfile(outDir,'TLS_checks_3controls.csv'));
writetable(B.Checks,fullfile(outDir,'TLS_checks_6controls.csv'));

% Matching check points: NOT controls in EITHER adjustment.
commonNames = intersect(A.Checks.Target,B.Checks.Target,'stable');
assert(~isempty(commonNames), 'No common independent check targets.');
[~,ia] = ismember(commonNames,A.Checks.Target);
[~,ib] = ismember(commonNames,B.Checks.Target);
T3 = A.Checks(ia,:); T6 = B.Checks(ib,:);
T = table(commonNames,T3.dX_mm,T3.dY_mm,T3.dZ_mm,T3.d3D_mm, ...
    T6.dX_mm,T6.dY_mm,T6.dZ_mm,T6.d3D_mm, ...
    'VariableNames',{'Target','dX_3_mm','dY_3_mm','dZ_3_mm','d3D_3_mm', ...
                     'dX_6_mm','dY_6_mm','dZ_6_mm','d3D_6_mm'});
T.Change6Minus3_mm = T.d3D_6_mm-T.d3D_3_mm;
writetable(T,fullfile(outDir,'TLS_checks_common_comparison.csv'));

% Report separate all-check summaries and directly comparable common checks.
labels = ["3 controls (all own checks)"; "6 controls (all own checks)"; ...
          "3 controls (common checks)"; "6 controls (common checks)"];
stats = [summaryStats(A.Checks);summaryStats(B.Checks); ...
         summaryStats(T3);summaryStats(T6)];
Summary = table(labels,stats(:,1),stats(:,2),stats(:,3),stats(:,4), ...
     stats(:,5),stats(:,6),stats(:,7),stats(:,8),stats(:,9), ...
     'VariableNames',{'Solution','N','RMSE_X_mm','RMSE_Y_mm', ...
        'RMSE_Z_mm','RMSE_3D_mm','Mean_dX_mm','Mean_dY_mm', ...
        'Mean_dZ_mm','Max_3D_mm'});
writetable(Summary,fullfile(outDir,'TLS_check_summary.csv'));

fprintf('\nCOMMON INDEPENDENT CHECK TARGETS (n=%d):\n',height(T));
disp(T);
fprintf('\nSUMMARY (mm):\n');
disp(Summary);
fprintf(['\nInterpret the last TWO rows when comparing the control schemes.\n' ...
  'Both rows then use exactly the same independently checked targets.\n']);

% One direct visual comparison of common-check 3D discrepancies.
f = figure('Name','Independent TLS check point comparison','Color','w');
b = bar(categorical(cellstr(commonNames)),[T.d3D_3_mm,T.d3D_6_mm]);
legend({'3 controls','6 controls'},'Location','best');
ylabel('3D check discrepancy (mm)');
xlabel('TLS check target');
title('Same independent targets in both adjustments');
grid on;
exportgraphics(f,fullfile(outDir,'TLS_common_checks_3D.png'),'Resolution',180);

% Coordinate-component discrepancies: signs help detect biases.
f2 = figure('Name','3D check displacement components','Color','w');
tiledlayout(3,1,'TileSpacing','compact');
fields = {'dX','dY','dZ'};
for j=1:3
    nexttile;
    switch j
        case 1, vals=[T.dX_3_mm,T.dX_6_mm];
        case 2, vals=[T.dY_3_mm,T.dY_6_mm];
        case 3, vals=[T.dZ_3_mm,T.dZ_6_mm];
    end
    bar(categorical(cellstr(commonNames)),vals); grid on;
    ylabel([fields{j} ' (mm)']); yline(0,'k-');
    if j==1,legend({'3 controls','6 controls'},'Location','best');end
end
sgtitle('Bundle minus TLS: independent check coordinate differences');
exportgraphics(f2,fullfile(outDir,'TLS_common_checks_components.png'),'Resolution',180);

fprintf('\nFiles saved in: %s\n',outDir);
end

function Out = extractChecks(matPath,tlsNames,tlsXYZ,label)
S = load(matPath);
assert(all(isfield(S,{'tie','tieRows','conRows','names','IDs'})), ...
    ['MAT file lacks fields expected from the Lab 5 script. ' ...
     'Need tie, tieRows, conRows, names, IDs.']);
assert(isfield(S.tie,'XYZ'),'Missing adjusted tie.XYZ.');
pointNames=upper(strtrim(string(S.names(:))));
pointIDs=double(S.IDs(:));
tieIDs=double(S.tieRows(:,1));
conIDs=double(S.conRows(:,1));
XYZadj=double(S.tie.XYZ);
assert(size(XYZadj,1)==numel(tieIDs) && size(XYZadj,2)==3, ...
    'Adjusted tie XYZ does not match tieRows ordering.');
assert(numel(pointNames)==numel(pointIDs), ...
    'names and IDs lengths disagree in MAT file.');
assert(numel(unique(pointIDs))==numel(pointIDs),'Duplicate bundle IDs.');
assert(isempty(intersect(tieIDs,conIDs)), ...
    'Control IDs unexpectedly occur in tieRows.');

% A name must be in TLS CSV, be part of the photogrammetric network, and
% be a TIE (not a control) in this adjustment.
[isInNetwork,ix] = ismember(tlsNames,pointNames);
isCheck=false(size(tlsNames));
idxInTie=zeros(size(tlsNames));
for i=find(isInNetwork(:))'
    [tf,k]=ismember(pointIDs(ix(i)),tieIDs);
    if tf && ~ismember(pointIDs(ix(i)),conIDs)
        isCheck(i)=true;
        idxInTie(i)=k;
    end
end
selected=find(isCheck);
assert(~isempty(selected),'No independent TLS checks in %s',label);
checkNames=tlsNames(selected);
adj=XYZadj(idxInTie(selected),:);
ref=tlsXYZ(selected,:);
delta=adj-ref;
d3=sqrt(sum(delta.^2,2));
Checks=table(checkNames,adj(:,1),adj(:,2),adj(:,3), ...
    ref(:,1),ref(:,2),ref(:,3), ...
    delta(:,1),delta(:,2),delta(:,3),d3, ...
    'VariableNames',{'Target','BundleX_mm','BundleY_mm','BundleZ_mm', ...
       'TLS_X_mm','TLS_Y_mm','TLS_Z_mm','dX_mm','dY_mm','dZ_mm','d3D_mm'});
Out.Checks=Checks;
Out.Label=label;
Out.ControlNames=pointNames(ismember(pointIDs,conIDs));
fprintf('%s: %d weighted controls, %d independent TLS checks\n', ...
    label,numel(conIDs),height(Checks));
end

function v=summaryStats(T)
err=[T.dX_mm,T.dY_mm,T.dZ_mm];
rmseComp=sqrt(mean(err.^2,1));
rmse3D=sqrt(mean(sum(err.^2,2)));
meanErr=mean(err,1);
v=[height(T),rmseComp,rmse3D,meanErr,max(T.d3D_mm)];
end
