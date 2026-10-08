% ENGO 531 Lab 2
close all;
clear;
clc;

%% Run adjustment

[eop,iop,con,tie,Cx,vr_img,vr_eop,vr_dis,vr_con] = ...
    full_matlab_bundle('engo531_lab2');


%% Network information

n_img = eop.num_eop;
n_tie = tie.num_tie;

n_img_obs = size(vr_img,1);
n_eop_obs = size(vr_eop,1);
n_dis_obs = size(vr_dis,1);


%% Question 7: n, u, r

n = n_img_obs + n_eop_obs + n_dis_obs;

u = 6*n_img + 3*n_tie;   % 6 EOP unknowns/image + XYZ/tie point

r = n-u;

% Check using redundancy numbers
r_check = sum(vr_img(:,2)) + ...
          sum(vr_eop(:,2)) + ...
          sum(vr_dis(:,2));

fprintf('n = %d\n',n)
fprintf('u = %d\n',u)
fprintf('r = %d\n',r)
fprintf('Sum of redundancy numbers = %.4f\n',r_check)


%% Calculate a-posteriori variance factor

% Observation standard deviations from CFG
sd_xy   = 0.7;       % pixels
sd_XYZ  = 0.0001;    % m = 0.1 mm
sd_wpk  = 0.0001;    % degrees
sd_dist = 0.010;     % m = 1 mm


% Image coordinate contribution
Q_img = sum((vr_img(:,1) ./ sd_xy).^2);


% EOP contribution
% vr_eop order:
% Xc, Yc, Zc, omega, phi, kappa

Q_eop_XYZ = sum((vr_eop(1:3,1) ./ sd_XYZ).^2);

Q_eop_wpk = sum((vr_eop(4:6,1) ./ sd_wpk).^2);


% Distance contribution
Q_dist = sum((vr_dis(:,1) ./ sd_dist).^2);


% Estimated a-posteriori variance factor
sigma0_sq = (Q_img + Q_eop_XYZ + Q_eop_wpk + Q_dist) / r;

% A-posteriori standard deviation of unit weight
sigma0 = sqrt(sigma0_sq);


fprintf('\nA-posteriori variance factor = %.4f\n',sigma0_sq)
fprintf('A-posteriori standard deviation = %.4f\n',sigma0)


%% Scale covariance matrix

% full_matlab_bundle returns an unscaled covariance matrix
Cx_scaled = sigma0_sq * Cx;


%% Question 8: Object point precision

n_eop_params = 6*n_img;
n_tie_params = 3*n_tie;


% Extract EOP covariance block
Cx_eops = Cx_scaled(1:n_eop_params,1:n_eop_params);


% Extract tie point covariance block
tie_start = n_eop_params + 1;

tie_end = n_eop_params + n_tie_params;

Cx_tie = Cx_scaled(tie_start:tie_end, ...
                   tie_start:tie_end);


%% Tie point standard deviations

tie_std_list = table( ...
    string(tie.pt(:)), ...
    zeros(n_tie,1), ...
    zeros(n_tie,1), ...
    zeros(n_tie,1), ...
    'VariableNames', ...
    {'Point','StdX_mm','StdY_mm','StdZ_mm'});


for i = 1:n_tie

    idx = 3*i-2:3*i;

    % 3x3 covariance matrix for tie point i
    Ci = Cx_tie(idx,idx);

    % sqrt(C) gives metres
    % multiply by 1000 to convert to mm
    tie_std_list.StdX_mm(i) = sqrt(Ci(1,1))*1000;
    tie_std_list.StdY_mm(i) = sqrt(Ci(2,2))*1000;
    tie_std_list.StdZ_mm(i) = sqrt(Ci(3,3))*1000;

end


% Combined 3D coordinate precision
tie_std_list.Std3D_mm = sqrt( ...
    tie_std_list.StdX_mm.^2 + ...
    tie_std_list.StdY_mm.^2 + ...
    tie_std_list.StdZ_mm.^2);


disp(tie_std_list)


%% Object point precision summary

fprintf('\nMean standard deviations:\n')

fprintf('X = %.3f mm\n', ...
    mean(tie_std_list.StdX_mm))

fprintf('Y = %.3f mm\n', ...
    mean(tie_std_list.StdY_mm))

fprintf('Z = %.3f mm\n', ...
    mean(tie_std_list.StdZ_mm))


% Best and worst determined points
[best_precision,best_idx] = min(tie_std_list.Std3D_mm);

[worst_precision,worst_idx] = max(tie_std_list.Std3D_mm);


fprintf('\nBest determined point:\n')
fprintf('Point %s, 3D precision = %.3f mm\n', ...
    tie_std_list.Point(best_idx), ...
    best_precision)


fprintf('\nWorst determined point:\n')
fprintf('Point %s, 3D precision = %.3f mm\n', ...
    tie_std_list.Point(worst_idx), ...
    worst_precision)


%% Error ellipsoids

figure;
hold on;
axis equal;
grid on;

S = 10;   % visual scale factor


for i = 1:n_tie

    idx = 3*i-2:3*i;

    Ci = Cx_tie(idx,idx);

    XYZi = tie.XYZ(i,:);

    calc_plot_ellipsoid(Ci,XYZi,S,20);

end


%% Plot camera positions

%{
plot3(eop.XYZc(:,1), ...
      eop.XYZc(:,2), ...
      eop.XYZc(:,3), ...
      'o', ...
      'MarkerSize',8, ...
      'LineWidth',1.5);


for i = 1:size(eop.XYZc,1)

    text(eop.XYZc(i,1), ...
         eop.XYZc(i,2), ...
         eop.XYZc(i,3), ...
         "Image " + string(eop.img(i)), ...
         'VerticalAlignment','bottom', ...
         'HorizontalAlignment','left');

end
%}

xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')

title('Object Point Error Ellipsoids Q9')

view(3)


%% Display EOP locations and orientations

EOP_table = table( ...
    eop.img(:), ...
    eop.XYZc(:,1), ...
    eop.XYZc(:,2), ...
    eop.XYZc(:,3), ...
    eop.wpk(:,1), ...
    eop.wpk(:,2), ...
    eop.wpk(:,3), ...
    'VariableNames', ...
    {'Image','Xc_m','Yc_m','Zc_m', ...
     'Omega_deg','Phi_deg','Kappa_deg'});

disp(EOP_table)