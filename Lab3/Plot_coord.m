% Photogrammetry Target Coordinates and EOPs
clear;
clc;
close all;

%% Import Target Coordinates
T = readmatrix('END 305 target coordinates transformed.txt');

ID = T(:,1);
X = T(:,2);
Y = T(:,3);
Z = T(:,4);

%% EOP Coordinates
% Columns: Image ID, X, Y, Z

EOP = [
    2   -280.7872   -4149.0255    152.2944;
    3   -137.8809   -3613.9113   1227.9180;
    4   1311.3245   -3712.2095    415.5383;
    5   2228.7159   -3592.4673   1293.5731;
    6   2228.7159   -2747.4673    253.5731;
    7  -1641.0098   -2594.2438   1587.8896;
    8   -162.8183   -2404.7969    683.2079
];


EOP_ID = EOP(:,1);
EOP_X = EOP(:,2);
EOP_Y = EOP(:,3);
EOP_Z = EOP(:,4);

%% Plot Coordinates
figure;
hold on;

% Plot target points
scatter3(X, Y, Z, 30, 'b', 'filled', ...
    'DisplayName', 'Target Points');

% Plot EOP locations
scatter3(EOP_X, EOP_Y, EOP_Z, 100, 'r', '^', 'filled', ...
    'DisplayName', 'EOP Locations');

%% Label Target Points
for i = 1:length(ID)
    text(X(i), Y(i), Z(i), num2str(ID(i)), ...
        'FontSize', 7);
end

%% Label EOP Locations
for i = 1:length(EOP_ID)
    text(EOP_X(i), EOP_Y(i), EOP_Z(i), ...
        sprintf(' Image %d', EOP_ID(i)), ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Color', 'r');
end

%% Plot Settings
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
title('Photogrammetry Target Points and EOP Locations');

legend('Location', 'best');
grid on;
axis equal;
view(3);
rotate3d on;
hold off;