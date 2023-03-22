T = 1;
dt = 2^(-10);

F = @(t, X) 0;
G = @(t, X) 1;
x0 = 0;
obj = sde(F, G, 'StartState', x0);    % dX = F(t,X)dt + G(t,X)dW

Nperiod = ceil(T/dt);
NTrials = 200;
for i = 1:length(powers)
    rng(24, 'twister')
    [x_smp,t_smp] = simByEuler(obj, Nperiod, 'DeltaTime', dt, 'NTrials', NTrials);
end

x_smp = squeeze(x_smp);


figure(1)
clf
hold on
% for i = 1:npow
    plot(t_smp, x_smp)
% end

title('Brownian Sample Paths', 'fontsize', 16)
xlabel('t')
ylabel('x(t)')

