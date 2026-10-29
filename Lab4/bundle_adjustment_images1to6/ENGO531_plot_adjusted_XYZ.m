%% ENGO 531 Lab 4 — Plot adjusted XYZ from a successful bundle adjustment
% Reads the *_results.mat saved by ENGO531_Lab4_select_images.m.
% Does not rerun or modify the adjustment.
clear; clc; close all;

[file,folder] = uigetfile('*_results.mat', ...
    'Choose the SUCCESSFUL six-image adjustment _results.mat');
if isequal(file,0), return; end
S = load(fullfile(folder,file));

assert(isfield(S,'tie') && isstruct(S.tie) && isfield(S.tie,'XYZ'), ...
    'The MAT file must contain tie.XYZ from the adjusted bundle.');
assert(isfield(S,'eop') && isstruct(S.eop) && isfield(S.eop,'XYZc'), ...
    'The MAT file must contain eop.XYZc.');
assert(isfield(S,'names') && isfield(S,'numericID'), ...
    'The MAT file must contain original point names and numericID.');

xyz = double(S.tie.XYZ);
cameras = double(S.eop.XYZc);
pointNames = string(S.names(:));
pointIDs = double(S.numericID(:));

assert(size(xyz,2)==3 && size(xyz,1)==numel(pointNames), ...
    'tie.XYZ must be one XYZ row per name.');
assert(size(cameras,2)==3, 'eop.XYZc must contain three columns.');
% Validate that adjusted coordinates follow the input numeric IDs.
if isfield(S.tie,'ID')
    assert(isequal(double(S.tie.ID(:)),pointIDs), ...
        'Tie-point IDs differ from input order; cannot safely label points.');
elseif isfield(S.tie,'id')
    assert(isequal(double(S.tie.id(:)),pointIDs), ...
        'Tie-point IDs differ from input order; cannot safely label points.');
end

% Image IDs follow rows of ext in the original script.
if isfield(S,'ext') && size(S.ext,1)==size(cameras,1)
    imageIDs = S.ext(:,1);
else
    imageIDs = (1:size(cameras,1))';
end

% Separate point groups (left and front in the successful six-image trial).
left = startsWith(pointNames,'L');
front = startsWith(pointNames,'F');
other = ~(left|front);

figure('Color','w','Name','ENGO 531 - Adjusted XYZ','Position',[100 100 1120 800]);
hold on;
if any(left)
    scatter3(xyz(left,1),xyz(left,2),xyz(left,3),60,[0.18 0.50 0.80],'filled', ...
        'DisplayName','Left-wall tie points');
end
if any(front)
    scatter3(xyz(front,1),xyz(front,2),xyz(front,3),60,[0.20 0.65 0.38],'filled', ...
        'DisplayName','Front-wall tie points');
end
if any(other)
    scatter3(xyz(other,1),xyz(other,2),xyz(other,3),60,[0.60 0.40 0.70],'filled', ...
        'DisplayName','Other tie points');
end
scatter3(cameras(:,1),cameras(:,2),cameras(:,3),125,[0.87 0.24 0.18], ...
    '^','filled','DisplayName','Camera stations');

for k=1:numel(pointNames)
    text(xyz(k,1),xyz(k,2),xyz(k,3),"  "+pointNames(k), ...
        'FontSize',9,'Color',[0.12 0.12 0.12]);
end
for k=1:numel(imageIDs)
    text(cameras(k,1),cameras(k,2),cameras(k,3), ...
        "  Image "+string(imageIDs(k)), 'FontSize',10,'FontWeight','bold', ...
        'Color',[0.70 0.15 0.12]);
end
% Connect camera positions in acquisition order.
plot3(cameras(:,1),cameras(:,2),cameras(:,3),'--', ...
    'Color',[0.87 0.24 0.18],'LineWidth',1,'HandleVisibility','off');

xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
title('Successful bundle adjustment — adjusted 3D tie points and EOP positions');
grid on; axis equal; view(35,25); rotate3d on;
legend('Location','bestoutside');

% Display numeric XYZ values, and optionally export figures/data.
pointTable = table(pointNames,pointIDs,xyz(:,1),xyz(:,2),xyz(:,3), ...
    'VariableNames',{'Target','BundleID','X_mm','Y_mm','Z_mm'});
cameraTable = table(imageIDs,cameras(:,1),cameras(:,2),cameras(:,3), ...
    'VariableNames',{'Image','Xc_mm','Yc_mm','Zc_mm'});
fprintf('\nADJUSTED TIE POINT COORDINATES [mm]:\n'); disp(pointTable);
fprintf('\nADJUSTED CAMERA COORDINATES [mm]:\n'); disp(cameraTable);

[~,base] = fileparts(file);
outputBase=fullfile(folder, [base '_3D_plot']);
exportgraphics(gcf,[outputBase '.png'],'Resolution',200);
writetable(pointTable,[outputBase '_points.csv']);
writetable(cameraTable,[outputBase '_cameras.csv']);
fprintf('\nSaved figure and CSVs to: %s*\n',outputBase);
