%SDE setup
%parameters from example 2 of https://arxiv.org/pdf/2208.10752.pdf

% b = -0.1;
sigma = 0.1;
% F = @(t,X) [X(2); -(X(1) +X(2) + 0.5*X(1)^3)];
F = @(t,X) [-2, -4; 5, -2]*X/2;
G = @(t,X) sigma * [0; 1];
% G = @(t,X) [0; 0.5*X(2)]/2;



%initial set
% x0 = [1; 1];
x0 = [0; 0.75];
obj = sde(F, G, 'StartState', x0);    % dX = F(t,X)dt + G(t,X)dW

dt = 1e-3;
T =5;

%Options
Nperiod = ceil(T/dt);
Antithetic = true;
NTrials = 1000;

[x_smp,t_smp] = simByEuler(obj, Nperiod, 'DeltaTime', dt, 'NTrials', NTrials,...
    'Antithetic', Antithetic);
% x_smp = squeeze(x_smp);

N_highlight = 4;
c = linspecer(4);

%% Plot


%unsafe set
theta_c = 5*pi/4; 
% Cu = [-1; -1];
% Ru = 0.3;

% theta_c = pi;
% theta_c = 3*pi/2;
Ru = 0.1;
% Cu = [-0.8; -0.8];
Cu = [-1;  -1];
% Cu = [-0.5;  -0.7];

figure(1)
clf
hold on
for i = 1:NTrials
%     plot3(t_smp, x_smp(:, 1, i), x_smp(:, 2, i));
%     plot(x_smp(:, 1, i), x_smp(:, 2, i));
if i > (NTrials - N_highlight)
    plot(x_smp(:, 1, i), x_smp(:, 2, i), 'color', c(NTrials-i+1, :), 'LineWidth', 2, 'HandleVisibility', 'off');
else
    plot(x_smp(:, 1, i), x_smp(:, 2, i), 'c', 'HandleVisibility', 'off');
end
end

% unsafe set plotting
theta_half_range = linspace(theta_c-pi/2, theta_c + pi/2, 200);
circ_half = [cos(theta_half_range); sin(theta_half_range)];
Xu = Cu + circ_half* Ru;
patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none', 'DisplayName', 'Unsafe Set')
       
%bounds from order-5 experiments/flow_test_std_distance_sde_point.m
% dist_mean = sqrt([2.59999995817759]);
% dist_eps = sqrt([1.51156001100412]);

dist_mean = sqrt([1.25055026768407]);
dist_cvar = sqrt([0.8064]);
dist_eps = sqrt([0.0234793503229038]);


x_dist_mean_align = dist_contour(100, Ru, dist_mean);
x_dist_cvar_align = dist_contour(100, Ru, dist_cvar);
x_dist_eps_align = dist_contour(100, Ru, dist_eps);

%transform the unsafe set
theta_cf = theta_c - 3*pi/2;
Rot_mat = [cos(theta_cf) -sin(theta_cf); sin(theta_cf) cos(theta_cf)];
x_dist_mean = Rot_mat*x_dist_mean_align + Cu;
x_dist_cvar = Rot_mat*x_dist_cvar_align + Cu;
x_dist_eps = Rot_mat*x_dist_eps_align + Cu;

    
plot(x_dist_mean(1, :), x_dist_mean(2, :), 'r--', 'DisplayName', 'Dist. Mean', 'LineWidth', 2)
plot(x_dist_cvar(1, :), x_dist_cvar(2, :), 'r:', 'DisplayName', 'Dist. CVAR 85%', 'LineWidth', 2)
plot(x_dist_eps(1, :), x_dist_eps(2, :), 'r', 'DisplayName', 'Dist. VP 85%', 'LineWidth', 2)


%TODO: plot the distance contour

scatter(x0(1), x0(2), 200, 'ko', 'HandleVisibility', 'off')
% titlestr = sprintf('dx = %0.2f x dt + %0.4f x dw', b, sigma);

xlabel('$x_1$', 'interpreter', 'latex')
ylabel('$x_2$', 'interpreter', 'latex')
title('Stochastic Distance Estimate', 'FontSize', 14)
legend('location', 'northwest')
pbaspect([diff(xlim), diff(ylim), 1])

function x_dist = dist_contour(Ntheta, R, c)
    %compute a contour at distance c away from the half-circle with N_theta
    %sample points


    theta_q1 = linspace(0, pi/2, Ntheta);
    theta_q2 = linspace(pi/2, pi, Ntheta);
    theta_q34 = linspace(pi, 2*pi, 2*Ntheta);

    %contour level
    

    x_right = [c*cos(theta_q1)+R; c*sin(theta_q1)];
    x_left = [c*cos(theta_q2)-R; c*sin(theta_q2)];
    x_bottom = [(c+R)*cos(theta_q34); (c+R)*sin(theta_q34)];
    x_dist = [x_right, x_left, x_bottom];
end
