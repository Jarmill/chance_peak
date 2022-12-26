%SDE setup
%parameters from example 1 of https://arxiv.org/pdf/2208.10752.pdf

a = 0;
b = 0;

T = 1;
F = @(t,X) (b-X)/(1-t);
%this is rational, we can deal with it through adding a new variable or
%using the Bugarin method https://arxiv.org/abs/1102.4954
G = @(t,X) 1;

x0 = a;
obj = sde(F, G, 'StartState', x0);    % dX = F(t,X)dt + G(t,X)dW

dt = 1e-3;


%Options
Nperiod = ceil(T/dt);
Antithetic = true;
NTrials = 100;

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