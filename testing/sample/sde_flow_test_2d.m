%SDE setup
%parameters from example 2 of https://arxiv.org/pdf/2208.10752.pdf

% b = -0.1;
sigma = 0.1;
F = @(t,X) [X(2); -(X(1) +X(2) + 0.5*X(1)^3)];
% F = @(t,X) [X(2); -(X(1) +X(2) - (1/3)*X(1)^3)];
% G = @(t,X) sigma * [0; X(2)];
G = @(t,X) sigma * [0; 1];


%unsafe set
theta_c = 5*pi/4; 
Cu = [-0.5; -0.75];
Ru = 0.5;

%initial set
x0 = [1; 1];
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
figure(1)
clf
hold on
for i = 1:NTrials
%     plot3(t_smp, x_smp(:, 1, i), x_smp(:, 2, i));
%     plot(x_smp(:, 1, i), x_smp(:, 2, i));
if i > (NTrials - N_highlight)
    plot(x_smp(:, 1, i), x_smp(:, 2, i), 'color', c(NTrials-i+1, :), 'LineWidth', 2);
else
    plot(x_smp(:, 1, i), x_smp(:, 2, i), 'c');
end
end

% unsafe set plotting
theta_half_range = linspace(theta_c-pi/2, theta_c + pi/2, 200);
circ_half = [cos(theta_half_range); sin(theta_half_range)];
Xu = Cu + circ_half* Ru;
patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none', 'DisplayName', 'Unsafe Set')
       
%TODO: plot the distance contour

scatter(x0(1), x0(2), 200, 'ko')
% titlestr = sprintf('dx = %0.2f x dt + %0.4f x dw', b, sigma);

xlabel('$x_1$', 'interpreter', 'latex')
ylabel('$x_2$', 'interpreter', 'latex')
title('Stochastic Flow System', 'FontSize', 14)
pbaspect([diff(xlim), diff(ylim), 1])


% plot(xlim, [1, 1]*[-0.874424669027551], 'r--', 'LineWidth', 3)
% plot(xlim, [1, 1]*-[1.16420596964602], 'r', 'LineWidth', 3)

% [0.874424669027551;1.16420596964602]
% 
% plot(xlim, [1, 1]*[-0.509913216820256], 'r--', 'LineWidth', 3)
% plot(xlim, [1, 1]*-[[0.950766054469061], 'r', 'LineWidth', 3)

% xlim([0, T]);