%SDE setup
%parameters from example 1 of https://arxiv.org/pdf/2208.10752.pdf
rng(25, 'twister');
b = -0.1;
sigma = sqrt(2)/2;
F = @(t,X) b * X;
G = @(t,X) sigma * X;

x0 = 0.25;
obj = sde(F, G, 'StartState', x0);    % dX = F(t,X)dt + G(t,X)dW

dt = 1e-3;
T = 1;

%Options
Nperiod = ceil(T/dt);
Antithetic = true;
% NTrials = 1000;
NTrials = 200;

[x_smp,t_smp] = simByEuler(obj, Nperiod, 'DeltaTime', dt, 'NTrials', NTrials,...
    'Antithetic', Antithetic);
x_smp = squeeze(x_smp);

%Plot
figure(1)
clf
hold on
plot(t_smp, x_smp);

titlestr = sprintf('dx = %0.2f x dt + %0.4f x dw', b, sigma);
title(titlestr, 'fontsize', 16)
xlabel('t')
ylabel('x(t)')

xlim([0, T]);