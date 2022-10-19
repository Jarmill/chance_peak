%SDE setup
%parameters from example 3 of https://arxiv.org/pdf/2208.10752.pdf

% b = -0.1;
sigma = 1;
F = @(t,X) [-2*X(1)+X(2)^2; -X(2)];
G = @(t,X) sigma*[1;1];

x0 = 0.5*[1; 1];
obj = sde(F, G, 'StartState', x0);    % dX = F(t,X)dt + G(t,X)dW

dt = 1e-3;
T =6;

%Options
Nperiod = ceil(T/dt);
Antithetic = true;
NTrials = 100;

[x_smp,t_smp] = simByEuler(obj, Nperiod, 'DeltaTime', dt, 'NTrials', NTrials,...
    'Antithetic', Antithetic);
% x_smp = squeeze(x_smp);

%% Plot
figure(1)
clf
hold on
for i = 1:NTrials
%     plot3(t_smp, x_smp(:, 1, i), x_smp(:, 2, i));
    plot(x_smp(:, 1, i), x_smp(:, 2, i));
end
scatter(x0(1), x0(2), 200, 'ko')
% titlestr = sprintf('dx = %0.2f x dt + %0.4f x dw', b, sigma);
titlestr='2d dynamics';
title(titlestr, 'fontsize', 16)
ylabel('x_1(t)')
ylabel('x_2(t)')

% xlim([0, T]);