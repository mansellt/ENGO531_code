%% ENGO 531 Lab 4 - Plot adjusted nine-image bundle with Image 5 datum
% Run after ENGO531_rebundle_datum_image5.m.
% Select ENGO531_joint_datum_image5_results.mat.
% Optionally select the original ENGO531_JOINT_1to9_point_map.csv to label targets.
clear; clc; close all;

[fn, folder] = uigetfile('*.mat','Select ENGO531_joint_datum_image5_results.mat');
if isequal(fn,0), return; end
S = load(fullfile(folder,fn));
assert(isfield(S,'tieResult') && isfield(S.tieResult,'XYZ'), ...
    'The selected file must contain tieResult.XYZ from the Image 5 adjustment.');
assert(isfield(S,'eop') && isfield(S.eop,'XYZc'), ...
    'The selected file must contain eop.XYZc.');
P=double(S.tieResult.XYZ);
C=double(S.eop.XYZc);
assert(size(P,2)==3 && size(C,2)==3,'XYZ must have three columns.');

% The course solver returns points in the original .tie row order.
if isfield(S,'T')
    pointID=double(S.T(:,1));
else
    pointID=(1:size(P,1))';
    warning('No input tie table in MAT; using row-number labels.');
end
assert(numel(pointID)==size(P,1),'Tie point IDs do not match adjusted XYZ rows.');

% Recover original camera numbering (bundle Image 1 = original Image 5).
if isfield(S,'origOrder')
    originalImages=double(S.origOrder(:));
else
    originalImages=[5;1;2;3;4;6;7;8;9];
end
assert(numel(originalImages)==size(C,1),'Camera count differs from image mapping.');

% Label targets using the point map produced by the earlier 9-image script.
labels="ID"+string(pointID);
mapFile=fullfile(folder,'ENGO531_JOINT_1to9_point_map.csv');
if ~isfile(mapFile)
    [mf,md]=uigetfile('*.csv', ...
        'Optional: select ENGO531_JOINT_1to9_point_map.csv (Cancel for numeric labels)',folder);
    if ~isequal(mf,0),mapFile=fullfile(md,mf);end
end
if isfile(mapFile)
    M=readtable(mapFile,'TextType','string');
    if all(ismember({'ID','Target'},M.Properties.VariableNames))
        [found,loc]=ismember(string(pointID),strtrim(string(M.ID)));
        labels(found)=string(M.Target(loc(found)));
        if ~all(found),warning('Some target IDs absent in point map; using numeric labels for those points.');end
    else
        warning('Point map lacks ID/Target columns. Using numeric labels.');
    end
end

colors=[0.16 0.43 0.84; 0.10 0.63 0.41; 0.91 0.44 0.15];
groups=["Left","Front","Right"];
mask={startsWith(labels,'L'),startsWith(labels,'F'),startsWith(labels,'R')};
others=~(mask{1}|mask{2}|mask{3});

fig=figure('Color','w','Name','Image 5 Datum - Adjusted 3D Network');
ax=axes(fig);hold(ax,'on');
for k=1:3
    idx=mask{k};
    if any(idx)
        scatter3(ax,P(idx,1),P(idx,2),P(idx,3),55,colors(k,:),'filled', ...
            'DisplayName',char(groups(k)));
    end
end
if any(others)
    scatter3(ax,P(others,1),P(others,2),P(others,3),50,[.55 .55 .55], ...
        'filled','DisplayName','Other targets');
end
scatter3(ax,C(:,1),C(:,2),C(:,3),105,'k','^','filled','DisplayName','Camera stations');
for i=1:size(P,1)
    text(ax,P(i,1),P(i,2),P(i,3)," "+labels(i),'FontSize',8,'Interpreter','none');
end
for i=1:size(C,1)
    text(ax,C(i,1),C(i,2),C(i,3),sprintf(' Image %d',originalImages(i)), ...
        'FontSize',9,'FontWeight','bold');
end
xlabel(ax,'X (mm)');ylabel(ax,'Y (mm)');zlabel(ax,'Z (mm)');
title(ax,'Adjusted 3D points and cameras - Original Image 5 datum');
grid(ax,'on');axis(ax,'equal');view(ax,35,25);rotate3d(fig,'on');
legend(ax,'Location','bestoutside');

fig2=figure('Color','w','Name','Image 5 Datum - Orthographic Views');
tiledlayout(fig2,1,3,'Padding','compact','TileSpacing','compact');
for j=1:3
    dims={[1 2],[1 3],[2 3]}; d=dims{j};
    ax=nexttile;hold(ax,'on');
    for k=1:3
        idx=mask{k};
        if any(idx)
            scatter(ax,P(idx,d(1)),P(idx,d(2)),35,colors(k,:),'filled');
        end
    end
    if any(others),scatter(ax,P(others,d(1)),P(others,d(2)),35,[.55 .55 .55],'filled');end
    scatter(ax,C(:,d(1)),C(:,d(2)),45,'k','^','filled');
    for i=1:size(P,1)
        text(ax,P(i,d(1)),P(i,d(2))," "+labels(i),'FontSize',7,'Interpreter','none');
    end
    xyzNames={'X','Y','Z'};
    xlabel(ax,[xyzNames{d(1)} ' (mm)']);ylabel(ax,[xyzNames{d(2)} ' (mm)']);
    title(ax,[xyzNames{d(1)} '-' xyzNames{d(2)} ' projection']);
    axis(ax,'equal');grid(ax,'on');
end

out=table(labels,pointID,P(:,1),P(:,2),P(:,3), ...
    'VariableNames',{'Target','PointID','X_mm','Y_mm','Z_mm'});
writetable(out,fullfile(folder,'ENGO531_image5_datum_adjusted_points.csv'));
exportgraphics(fig,fullfile(folder,'ENGO531_image5_datum_3D.png'),'Resolution',220);
exportgraphics(fig2,fullfile(folder,'ENGO531_image5_datum_views.png'),'Resolution',220);
fprintf('Plotted %d adjusted points and %d cameras.\n',size(P,1),size(C,1));
fprintf('Saved plots and coordinate CSV in %s\n',folder);
