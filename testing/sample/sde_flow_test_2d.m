%SDE setup
%parameters from example 2 of https://arxiv.org/pdf/2208.10752.pdf

% b = -0.1;
sigma = 0.1;
F = @(t,X) [X(2); -(X(1) +X(2) + 0.5*X(1)^3)];
% F = @(t,X) [X(2); -(X(1) +X(2) - (1/3)*X(1)^3)];
% G = @(t,X) sigma * [0; X(2)];
G = @(t,X) sigma * [0; 1];

x0 = [1; 1];
obj = sde(F, G, 'StartState', x0);    % dX = F(t,X)dt + G(t,X)dW

dt = 1e-3;
T =5;

%Options
Nperiod = ceil(T/dt);
Antithetic = true;
% NTrials = 10000;
NTrials = 500;

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

if i > (NTrials - N_highlight)
    plot(x_smp(:, 1, i), x_smp(:, 2, i), 'color', c(NTrials-i+1, :), 'LineWidth', 2);
else
    plot(x_smp(:, 1, i), x_smp(:, 2, i), 'c');
end
end
scatter(x0(1), x0(2), 200, 'ko')

xlabel('$x_1$', 'interpreter', 'latex')
ylabel('$x_2$', 'interpreter', 'latex')
title('Stochastic Flow System', 'FontSize', 14)

plot(xlim, [1, 1]*[-0.874424669027551], 'r--', 'LineWidth', 3)
plot(xlim, [1, 1]*-[1.16420596964602], 'r', 'LineWidth', 3)

% [0.874424669027551;1.16420596964602]
% 
% plot(xlim, [1, 1]*[-0.509913216820256], 'r--', 'LineWidth', 3)
% plot(xlim, [1, 1]*-[[0.950766054469061], 'r', 'LineWidth', 3)

% xlim([0, T]);