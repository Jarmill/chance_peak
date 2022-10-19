%SDE setup
%parameters from example 5 of https://arxiv.org/pdf/2208.10752.pdf

% b = -0.1;
sigma = 1;

omega=1;
a0 = 1;
alpha0=1;
T0=1;
eta = 0.5;
F = @(t,X) [X(2); -omega^2*X(1) - a0*X(2)];
G = @(t,X) [0;eta*(-alpha0*X(2) + T0)];

x0 = 0.25*[1; 1];
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
xlabel('x_1(t)')
ylabel('x_2(t)')

% xlim([0, T]);