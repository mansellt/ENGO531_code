% Photogrammetry Target Coordinates
clear;
clc;
close all;

% Import coordinates
T = readmatrix('END 305 target coordinates transformed.txt');

% Extract columns
ID = T(:,1);
X = T(:,2);
Y = T(:,3);
Z = T(:,4);

% 3D plot
figure;
scatter3(X, Y, Z, 40, Z, 'filled');
hold on;

% Label points
for i = 1:length(ID)
    text(X(i), Y(i), Z(i), num2str(ID(i)), ...
        'FontSize', 8);
end

xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
title('Photogrammetry Target Coordinates');

colorbar;
grid on;
axis equal;
view(3);