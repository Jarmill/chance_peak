T = 1;
powers = 2:9;
N = 2.^powers;

% Nperiod = T/N;
npow = length(powers);
x_smp = cell(npow, 1);
t_smp  = cell(npow, 1);

F = @(t, X) 0;
G = @(t, X) 1;
x0 = 0;
obj = sde(F, G, 'StartState', x0);    % dX = F(t,X)dt + G(t,X)dW


for i = 1:length(powers)
    rng(24, 'twister')
    [x_smp{i},t_smp{i}] = simByEuler(obj, N(i), 'DeltaTime', 1/N(i), 'NTrials', 1);
end


figure(1)
clf
hold on
for i = 1:npow
    plot(t_smp{i}, x_smp{i})
end

title('Improved Brownian Discretization', 'fontsize', 16)
xlabel('t')
ylabel('x(t)')

