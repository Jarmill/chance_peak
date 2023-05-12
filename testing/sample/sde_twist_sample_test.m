%SDE setup
%parameters from example 2 of https://arxiv.org/pdf/2208.10752.pdf

% b = -0.1;
sigma = 0.1;
A_true = [-1 1 1; -1 0 -1; 0 1 -2];
B_true = [-1 0 -1;
          0 1 1;
          1 1 0]/2;

F = @(t,x) A_true*x - B_true*(4*x.^3 - 3*x);
% F = @(t,X) [X(2); -(X(1) +X(2) - (1/3)*X(1)^3)];
% G = @(t,X) sigma * [0; X(2)];
% G = @(t,x) sigma * [0; 0; 1];
G = @(t,x) sigma * [0; 0;1];

x0 = [-0.5; 0; 0];
obj = sde(F, G, 'StartState', x0);    % dX = F(t,X)dt + G(t,X)dW

dt = 1e-3;
T =5;

%Options
Nperiod = ceil(T/dt);
Antithetic = true;
% NTrials = 500;
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
    plot3(x_smp(:, 1, i), x_smp(:, 2, i),x_smp(:, 3, i), 'color', c(NTrials-i+1, :), 'LineWidth', 2);
else
    plot3(x_smp(:, 1, i), x_smp(:, 2, i),x_smp(:, 3, i), 'c');
end
end
scatter3(x0(1), x0(2), x0(3), 200, 'ko')
% titlestr = sprintf('dx = %0.2f x dt + %0.4f x dw', b, sigma);

xlabel('$x_1$', 'interpreter', 'latex')
ylabel('$x_2$', 'interpreter', 'latex')
zlabel('$x_3$', 'interpreter', 'latex')
title('Stochastic Twist System', 'FontSize', 14)
view(3)


%% plot the patches
% chance_p = [0.820253508023993;0.975507232859090];
%degree 6 bounds
chance_p = [0.820253508023993; %mean
    0.973259607789189; %CVAR
    1.32019648246837]; %VP
xl = xlim;
yl = ylim;
patch(xl([1,1,2,2,1]), yl([1,2,2,1,1]), chance_p(1)*ones(1, 5), 'r', 'EdgeColor', 'None', 'FaceAlpha', 0.5)
patch(xl([1,1,2,2,1]), yl([1,2,2,1,1]), chance_p(2)*ones(1, 5), 'k', 'EdgeColor', 'None', 'FaceAlpha', 0.5)
patch(xl([1,1,2,2,1]), yl([1,2,2,1,1]), chance_p(3)*ones(1, 5), 'r', 'EdgeColor', 'None')


% plot(xlim, [1, 1]*[-0.874424669027551], 'r--', 'LineWidth', 3)
% plot(xlim, [1, 1]*-[1.16420596964602], 'r', 'LineWidth', 3)

% [0.874424669027551;1.16420596964602]
% 
% plot(xlim, [1, 1]*[-0.509913216820256], 'r--', 'LineWidth', 3)
% plot(xlim, [1, 1]*-[[0.950766054469061], 'r', 'LineWidth', 3)

% xlim([0, T]);