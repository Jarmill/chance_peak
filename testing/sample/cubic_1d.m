%SDE setup
%parameters from example 1 of https://arxiv.org/pdf/2208.10752.pdf
rng(25, 'twister');
% b = -0.1;
sigma = 0.2;
F = @(t,X) X - X.^3;
% G = @(t,X) sigma * X;
G = @(t,X) sigma*ones(size(X));

% x0 = 0.5;

N0 = 41;
x0 = linspace(-0.75, 0.75, N0)';
% x0 = 0;
obj = sde(F, G, 'StartState', x0);    % dX = F(t,X)dt + G(t,X)dW

dt = 1e-3;
T = 1;
% T = 5;
% T = 20;
% T = 50;

%Options
Nperiod = ceil(T/dt);
Antithetic = true;
% NTrials = 1000;
% NTrials = 500;
% NTrials = 200;
% NTrials = 5;
NTrials = 100;

[x_smp,t_smp] = simByEuler(obj, Nperiod, 'DeltaTime', dt, 'NTrials', NTrials,...
    'Antithetic', Antithetic);
x_smp = squeeze(x_smp);

%Plot
figure(1)
clf
hold on
if N0 == 1
plot(t_smp, x_smp);
else
    for i = 1:N0
        plot(t_smp, squeeze(x_smp(:, i, :)));
    end
end

titlestr = sprintf('dx = (x-x^3) dt + %0.2f x dw', sigma);
title(titlestr, 'fontsize', 16)
xlabel('t')
ylabel('x(t)')

xlim([0, T]);