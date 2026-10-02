%ENGO 531 Lab1
close all;
clear all;
clc;

%% Run Adjustment
% Original adjustment sd_xy = 0.3 pix
[eop1,iop1,con1,tie1,Cx1,vr_img1,vr_eop1,vr_dis1,vr_con1]=full_matlab_bundle('engo531_lab1');

% Improved measurement simulation: sd_xy = 0.03 px
[eop2,iop2,con2,tie2,Cx2,vr_img2,vr_eop2,vr_dis2,vr_con2] = full_matlab_bundle('engo531_lab1_2');

% improved sd_xy = 0.03 pix
n_img_xy = size(vr_img1, 1); %number of image measurements
n_con = size(vr_con1, 1); %XYZ for each control

%check residuals of image coordinates
%figure;
%scatter(1:n_img_xy, vr_img(:,1), "filled", "o" )
%figure;
%scatter(1:n_con, vr_con(:,1), "filled", "o" )

%        Cx matrix Structure
%4 images with 6 EOPS each: 24
%8 Tiepoints with XYZ = 24
%4 Control with XYZ = 12
% Cx (60 x 60)
% XYZ in mm and omega phi kappa in decimal deg

%get 3x3 covariance-variance matrix for each control/tie point

Cx_tie_ctr = Cx1(25:60, 25:60); %just take submatrix of tie and control points

point_list_str = [tie1.pt;    %list of points
                  con1.pt];

XYZ_tie_ctr = [tie1.XYZ; %XYZ coordinatres of tie and control
               con1.XYZ];

n_points = length(XYZ_tie_ctr); %number of points


%% 3. Plot and Analyze results
figure;
hold on;
axis equal;
grid on;

S=10; %scale factor

tie_ctr_std_list = table( ...
    point_list_str(:), zeros(n_points,1), zeros(n_points,1), zeros(n_points,1), ...
    'VariableNames', {'Point','StdX_mm','StdY_mm','StdZ_mm'});

for i = 1:n_points

    %118.555234944390 %From simulated adjustment used to scale covariance

    idx = 3*i-2:3*i; % get the 3x3 covariance matrix for a point
    Ci = Cx_tie_ctr(idx,idx);
    XYZi = XYZ_tie_ctr(i,:); %point coordinates [tie; ctr]
    calc_plot_ellipsoid(Ci, XYZi, S, 20); % Plot ellipsoid
    text(XYZi(1), XYZi(2), XYZi(3), point_list_str(i),'VerticalAlignment','bottom','HorizontalAlignment','left');

    tie_ctr_std_list.StdX_mm(i) = sqrt(Ci(1,1));
    tie_ctr_std_list.StdY_mm(i) = sqrt(Ci(2,2));
    tie_ctr_std_list.StdZ_mm(i) = sqrt(Ci(3,3));


end

%{
%% Plot the Camera perspective
plot3(eop1.XYZc(:,1), eop1.XYZc(:,2), eop1.XYZc(:,3), 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
for i = 1:size(eop1.XYZc,1)
    text(eop1.XYZc(i,1), eop1.XYZc(i,2), eop1.XYZc(i,3), string(eop1.img(i)), ...
         'VerticalAlignment','bottom', 'HorizontalAlignment','left');
end
%}


xlabel('X mm');
ylabel('Y mm');
zlabel('Z mm');
title('Point 3D Error Ellipsoids');
view(3)

tie_ctr_std_list

%standard deviation for EOPS
Cx_eops = Cx1(1:24, 1:24);
n_img = length(eop1.img);  %number of images
eop_std_list = table( ...
    eop1.img(:), zeros(n_img,1), zeros(n_img,1), zeros(n_img,1),zeros(n_img,1), ... 
    zeros(n_img,1),zeros(n_img,1) , ...
    'VariableNames', {'img_ID','StdX_mm','StdY_mm','StdZ_mm', ...
    'StdOmega_DD','StdPhi_DD','StdKappa_DD'});

%build std list for EOPs
for i = 1:n_img
    idx = 6*i-5:6*i; %get 6x6 matrix
    Ci = Cx_eops(idx,idx);
    eop_std_list.StdX_mm(i) = sqrt(Ci(1,1));
    eop_std_list.StdY_mm(i) = sqrt(Ci(2,2));
    eop_std_list.StdZ_mm(i) = sqrt(Ci(3,3));
    eop_std_list.StdOmega_DD(i) = sqrt(Ci(4,4));
    eop_std_list.StdPhi_DD(i) = sqrt(Ci(5,5));
    eop_std_list.StdKappa_DD(i) = sqrt(Ci(6,6));
end 

eop_std_list

%%  /////////////////     Questions       ////////////////////////////////////////////////////
%% a. What are n, u and r for this network? r can be computed two ways. Explain each.

    %r(dof) can be computed by summing the redundancy numbers of each observation
r = sum(vr_img1(:,2)) + sum(vr_con1(:, 2));

    %or by n-u
n = n_img_xy + n_con;
u = tie1.num_tie*3 + con1.num_con*3 + eop1.num_eop*6; % XYZ(control and tie points) XYZwphikappa

r_check = n-u;

%% b. Focusing on the image point observations only, are there any outliers? Explain.
point_image = [ %This is the vector of the point number and image number
     1     3
     2     3
     3     3
     4     3
     5     3
     9     3
    12     3
    13     3
    15     3
    20     3
     2    10
     3    10
     4    10
     5    10
    12    10
    15    10
    20    10
     1    12
     2    12
     3    12
     4    12
     9    12
    15    12
    17    12
    20    12
    21    12
     1    13
     2    13
     3    13
     4    13
     5    13
     9    13
    13    13
    15    13
    17    13
    21    13
];

w_img = vr_img1(:,1) ./ (0.3 .* sqrt(vr_img1(:,2)));

crit = 1.96;

outlier_idx = find(abs(w_img) > crit);
w_img_outlier = w_img(outlier_idx);
w_xy = reshape(w_img,2,[])';

outlier_table = table(point_image(:,1), point_image(:,2), w_xy(:,1), w_xy(:,2), ...
    'VariableNames', {'Point','Image','wx','wy'});

outlier_table.OutlierX = abs(outlier_table.wx) > 1.96;
outlier_table.OutlierY = abs(outlier_table.wy) > 1.96;

%% c. Remembering that there are two sets of observations in this network, compute the estimated
%variance factor (EVF). Analyze the EVF. What conclusion can be drawn from it?

sd_xy1 = 0.3; %px for first adjustment
sd_xy2 = 0.03; %px for second adjustment
sd_XYZ = 1;  %mm

comb_var1 = (vr_img1(:,1)'*(sd_xy1^-2)*vr_img1(:,1) + vr_con1(:,1)'*(sd_XYZ^-2)*vr_con1(:,1))/r
comb_var2 = (vr_img2(:,1)'*(sd_xy2^-2)*vr_img2(:,1) + vr_con2(:,1)'*(sd_XYZ^-2)*vr_con2(:,1))/r
%EVH = 1: assumed observation precision is consistent with residuals
%EVH > 1: residuals are larger than the expected a-priori precision
%EVH < 1: residuals are smaller than the expected a-priori precision

        
%% d. Focusing on the tie points, which point is the most strongly determined? Which is the most
%weakly determined? Justify your response in both cases. Comment on the overall tie point
%precision.


% using trace of sub variance covariance matrix for each point
tie_std_list = tie_ctr_std_list(1:8,:);

sum_var = tie_std_list.StdX_mm.^2+tie_std_list.StdY_mm.^2+ tie_std_list.StdZ_mm.^2;

[min_sum_var min_idx]= min(sum_var(:,1));
[max_sum_var max_idx]= max(sum_var(:,1));

max_tie = tie_std_list(max_idx,"Point");
min_tie = tie_std_list(min_idx,"Point");
%21 is weakest, 4 is strongest

% Also use error ellipse semi major and semi minor axis

%% e. Focusing again on the image points, which one is the least controlled? Which is the most
%controlled? Explain why for each case and the ramification(s).

    % least controlled = lowest redundancy number
    %  most controlled = highest redundancy number
    % Look at redundancy for each coordinate
    % look at point in image: So plot camera coordinates as well

r_xy = reshape(vr_img1(:,2), 2, [])';
% Find overall smallest and largets redundancy
[max_r, linear_idx_max] = max(r_xy(:));
[min_r, linear_idx_min] = min(r_xy(:)); 

[row_max, col_max] = ind2sub(size(r_xy), linear_idx_max); % Convert linear index into row and column
[row_min, col_min] = ind2sub(size(r_xy), linear_idx_min);

max_img_point = point_image(row_max,:); % Get point number and image number
min_img_point = point_image(row_min,:);

% Determine whether it was x or y
coord = ["x","y"];

max_coord = coord(col_max);
min_coord = coord(col_min);

% Display
fprintf('Most controlled: Point %d in Image %d, %s-coordinate, r = %.4f\n', ...
    max_img_point(1), max_img_point(2), max_coord, max_r);

fprintf('Least controlled: Point %d in Image %d, %s-coordinate, r = %.4f\n', ...
    min_img_point(1), min_img_point(2), min_coord, min_r);

%Max: Point 15 in image 1

%% f. Perform a simulation using the same input data. Re-run the adjustment with a smaller
%standard deviation for the image point observations. Use 0.03 pix to simulate the availability
%of a superior measurement algorithm. What is/are the effect(s)? What is the main conclusion
%that can be drawn from this simulation?

                                % Original adjustment outputs:
                                % Cx1, tie1, con1
                                % Improved adjustment outputs:
                                % Cx2, tie2, con2

% Tie point covariance blocks
Cx_tie1 = Cx1(25:48, 25:48);
Cx_tie2 = Cx2(25:48, 25:48);

n_tie = size(tie1.XYZ,1);

tie_compare = table( ...
    string(tie1.pt(:)), ...
    zeros(n_tie,1), zeros(n_tie,1), zeros(n_tie,1), ...
    zeros(n_tie,1), zeros(n_tie,1), zeros(n_tie,1), ...
    'VariableNames', {'Point', ...
    'StdX_old','StdY_old','StdZ_old', ...
    'StdX_new','StdY_new','StdZ_new'});

for i = 1:n_tie

    idx = 3*i-2:3*i;

    C1 = Cx_tie1(idx,idx);
    C2 = Cx_tie2(idx,idx);

    tie_compare.StdX_old(i) = sqrt(C1(1,1));
    tie_compare.StdY_old(i) = sqrt(C1(2,2));
    tie_compare.StdZ_old(i) = sqrt(C1(3,3));

    tie_compare.StdX_new(i) = sqrt(C2(1,1));
    tie_compare.StdY_new(i) = sqrt(C2(2,2));
    tie_compare.StdZ_new(i) = sqrt(C2(3,3));
end


% Control point comparison
% Control points are after EOPs and tie points
Cx_con1 = Cx1(49:60, 49:60);
Cx_con2 = Cx2(49:60, 49:60);

n_con = size(con1.XYZ,1);

con_compare = table( ...
    string(con1.pt(:)), ...
    zeros(n_con,1), zeros(n_con,1), zeros(n_con,1), ...
    zeros(n_con,1), zeros(n_con,1), zeros(n_con,1), ...
    'VariableNames', {'Point', ...
    'StdX_old','StdY_old','StdZ_old', ...
    'StdX_new','StdY_new','StdZ_new'});

for i = 1:n_con

    idx = 3*i-2:3*i;

    C1 = Cx_con1(idx,idx);
    C2 = Cx_con2(idx,idx);

    con_compare.StdX_old(i) = sqrt(C1(1,1));
    con_compare.StdY_old(i) = sqrt(C1(2,2));
    con_compare.StdZ_old(i) = sqrt(C1(3,3));

    con_compare.StdX_new(i) = sqrt(C2(1,1));
    con_compare.StdY_new(i) = sqrt(C2(2,2));
    con_compare.StdZ_new(i) = sqrt(C2(3,3));
end


%% Tie point coordinate changes
tie_coord_change = table( ...
    string(tie1.pt(:)), ...
    tie2.XYZ(:,1) - tie1.XYZ(:,1), ...
    tie2.XYZ(:,2) - tie1.XYZ(:,2), ...
    tie2.XYZ(:,3) - tie1.XYZ(:,3), ...
    'VariableNames', {'Point','dX_mm','dY_mm','dZ_mm'});

tie_coord_change.d3D_mm = sqrt( ...
    tie_coord_change.dX_mm.^2 + ...
    tie_coord_change.dY_mm.^2 + ...
    tie_coord_change.dZ_mm.^2);

disp(tie_coord_change)


%% Plot coordinate changes - Tie Points

figure;

bar(categorical(tie_coord_change.Point), ...
    [tie_coord_change.dX_mm, ...
     tie_coord_change.dY_mm, ...
     tie_coord_change.dZ_mm]);

xlabel('Tie Point');
ylabel('Coordinate Change (mm)');
title('Tie Point Coordinate Changes: 0.3 px vs 0.03 px');
legend('\DeltaX','\DeltaY','\DeltaZ','Location','best');
grid on;
yline(0,'k-');

%% Control point coordinate changes

con_coord_change = table( ...
    string(con1.pt(:)), ...
    con2.XYZ(:,1) - con1.XYZ(:,1), ...
    con2.XYZ(:,2) - con1.XYZ(:,2), ...
    con2.XYZ(:,3) - con1.XYZ(:,3), ...
    'VariableNames', {'Point','dX_mm','dY_mm','dZ_mm'});

con_coord_change.d3D_mm = sqrt( ...
    con_coord_change.dX_mm.^2 + ...
    con_coord_change.dY_mm.^2 + ...
    con_coord_change.dZ_mm.^2);

%% Plot coordinate changes - Control Points

figure;

bar(categorical(con_coord_change.Point), ...
    [con_coord_change.dX_mm, ...
     con_coord_change.dY_mm, ...
     con_coord_change.dZ_mm]);

xlabel('Control Point');
ylabel('Coordinate Change (mm)');
title('Control Point Coordinate Changes: 0.3 px vs 0.03 px');
legend('\DeltaX','\DeltaY','\DeltaZ','Location','best');
grid on;
yline(0,'k-');

%% Total 3D coordinate change - Tie Points

figure;

bar(categorical(tie_coord_change.Point), ...
    tie_coord_change.d3D_mm);

xlabel('Tie Point');
ylabel('3D Coordinate Change (mm)');
title('Magnitude of Tie Point Coordinate Changes');
grid on;


%% Total 3D coordinate change - Control Points

figure;

bar(categorical(con_coord_change.Point), ...
    con_coord_change.d3D_mm);

xlabel('Control Point');
ylabel('3D Coordinate Change (mm)');
title('Magnitude of Control Point Coordinate Changes');
grid on;