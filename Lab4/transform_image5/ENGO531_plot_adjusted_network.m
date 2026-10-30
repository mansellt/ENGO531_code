%% ENGO 531 Lab 4 - Plot the ADJUSTED 9-image photogrammetric network
% Uses the converged _results.mat and the matching _point_map.csv.
% Does not rerun the bundle adjustment or require earlier scripts.
clear; clc; close all;

%% Options
showPointLabels = true;
showCameraLabels = true;
showCameras = true;
showCameraToObjectLines = false;  % set true for a rough network overview
saveFigures = true;
exportAdjustedCSV = true;

%% Select files
[fn, folder] = uigetfile('*.mat', 'Select ENGO531_JOINT_1to9_results.mat');
if isequal(fn,0), return; end
matPath = fullfile(folder, fn);
S = load(matPath);
assert(isfield(S,'tie') && isfield(S.tie,'XYZ'), ...
    'Results MAT is missing tie.XYZ (adjusted object coordinates).');
assert(isfield(S,'eop') && isfield(S.eop,'XYZc'), ...
    'Results MAT is missing eop.XYZc (adjusted camera coordinates).');
P = double(S.tie.XYZ);
C = double(S.eop.XYZc);
assert(size(P,2)==3 && size(C,2)==3, 'XYZ data must have three columns.');

mapPath = fullfile(folder,'ENGO531_JOINT_1to9_point_map.csv');
if ~isfile(mapPath)
    [mf, md] = uigetfile('*.csv', 'Select matching ENGO531_JOINT_1to9_point_map.csv',folder);
    if isequal(mf,0), return; end
    mapPath = fullfile(md,mf);
end
M = readtable(mapPath, 'TextType','string');
assert(all(ismember({'Target','ID'},M.Properties.VariableNames)), ...
    'Point map must contain columns Target and ID.');
assert(height(M)==size(P,1), ...
    'Map has %d targets, but adjusted result has %d points.',height(M),size(P,1));

% The course output tie.XYZ is in its tie-file row order. Check point map
% against the original bundle .tie file when it is present.
originalTiePath = fullfile(folder,'ENGO531_JOINT_1to9.tie');
if isfile(originalTiePath)
    T0 = readmatrix(originalTiePath,'FileType','text');
    assert(size(T0,1)==size(P,1), 'Original .tie count differs from result.');
    assert(all(T0(:,1)==M.ID), ...
        'Point-map ID order differs from original .tie order. Cannot safely label results.');
else
    warning('Original .tie not found: assuming point-map rows match tie.XYZ row ordering.');
end

labels = string(M.Target);
imgID = (1:size(C,1))';
if isfield(S,'ext') && isnumeric(S.ext) && size(S.ext,1)==size(C,1)
    imgID = double(S.ext(:,1));
end

groups = strings(numel(labels),1);
groups(startsWith(labels,'L','IgnoreCase',true)) = "Left";
groups(startsWith(labels,'F','IgnoreCase',true)) = "Front";
groups(startsWith(labels,'R','IgnoreCase',true)) = "Right";
assert(all(groups~=""), 'Some point names do not begin with L, F, or R.');
colors = [0.15 0.42 0.82; 0.10 0.63 0.40; 0.88 0.41 0.14];
groupNames = ["Left","Front","Right"];

%% 3D view
fig3 = figure('Color','w','Name','Adjusted 3D Network','Position',[100 100 1100 760]);
ax=axes(fig3); hold(ax,'on');
for k=1:3
    sel=groups==groupNames(k);
    scatter3(ax,P(sel,1),P(sel,2),P(sel,3),55,colors(k,:), ...
        'filled','DisplayName',groupNames(k)+" targets");
    if showPointLabels
        ii=find(sel);
        for j=ii(:)'
            text(ax,P(j,1),P(j,2),P(j,3),"  "+labels(j), ...
                'Color',colors(k,:),'FontSize',9,'Interpreter','none');
        end
    end
end
if showCameras
    scatter3(ax,C(:,1),C(:,2),C(:,3),100,'k','^','filled','DisplayName','Camera stations');
    if showCameraLabels
        for j=1:size(C,1)
            text(ax,C(j,1),C(j,2),C(j,3),sprintf('  Image %d',imgID(j)), ...
                'Color','k','FontWeight','bold','FontSize',9);
        end
    end
    if showCameraToObjectLines
        mid=mean(P,1);
        for j=1:size(C,1)
            plot3(ax,[C(j,1) mid(1)],[C(j,2) mid(2)],[C(j,3) mid(3)], ...
                ':','Color',[0.75 0.75 0.75],'HandleVisibility','off');
        end
    end
end
xlabel(ax,'X (mm)'); ylabel(ax,'Y (mm)'); zlabel(ax,'Z (mm)');
title(ax,'Adjusted object points and camera stations (9-image bundle)');
grid(ax,'on'); axis(ax,'equal'); view(ax,35,25); rotate3d(fig3,'on');
legend(ax,'Location','bestoutside');

%% Orthographic views: useful for checking individual faces
fig2 = figure('Color','w','Name','Adjusted Plan and Elevations','Position',[140 140 1200 760]);
tiledlayout(fig2,1,3,'Padding','compact','TileSpacing','compact');
views = {[1 2],[1 3],[2 3]};
xl = {'X (mm)','X (mm)','Y (mm)'};
yl = {'Y (mm)','Z (mm)','Z (mm)'};
ttl = {'Plan: X-Y','Front projection: X-Z','Side projection: Y-Z'};
for v=1:3
    ax = nexttile; hold(ax,'on'); dims=views{v};
    for k=1:3
        sel=groups==groupNames(k);
        scatter(ax,P(sel,dims(1)),P(sel,dims(2)),32,colors(k,:),'filled', ...
            'DisplayName',groupNames(k));
        if showPointLabels
            ii=find(sel);
            for j=ii(:)'
                text(ax,P(j,dims(1)),P(j,dims(2))," "+labels(j), ...
                    'Color',colors(k,:),'FontSize',8,'Interpreter','none');
            end
        end
    end
    if showCameras
        scatter(ax,C(:,dims(1)),C(:,dims(2)),48,'k','^','filled','DisplayName','Cameras');
    end
    xlabel(ax,xl{v}); ylabel(ax,yl{v}); title(ax,ttl{v});
    axis(ax,'equal'); grid(ax,'on');
end

%% Export adjusted points
adjustedPoints=table(labels,M.ID,groups,P(:,1),P(:,2),P(:,3), ...
    'VariableNames',{'Target','ID','Face','X_mm','Y_mm','Z_mm'});
disp('Adjusted object-point coordinates [mm]:');
disp(adjustedPoints);
if exportAdjustedCSV
    csvOut=fullfile(folder,'ENGO531_JOINT_1to9_ADJUSTED_XYZ.csv');
    writetable(adjustedPoints,csvOut);
    fprintf('Saved adjusted XYZ table: %s\n',csvOut);
end
if saveFigures
    png3=fullfile(folder,'ENGO531_adjusted_3D_network.png');
    png2=fullfile(folder,'ENGO531_adjusted_orthographic_views.png');
    exportgraphics(fig3,png3,'Resolution',200);
    exportgraphics(fig2,png2,'Resolution',200);
    fprintf('Saved figures:\n  %s\n  %s\n',png3,png2);
end
