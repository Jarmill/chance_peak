rng(40, 'twister')
dt = 1e-3;
T =5;


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
p = @(x) squeeze(x(:, 3, :));

x0 = [-0.5; 0; 0];
obj = sde(F, G, 'StartState', x0);    % dX = F(t,X)dt + G(t,X)dW



Nperiod = ceil(T/dt);
Antithetic = true;
NTrials = 50000;
Nblock = 1000;
p_smp = zeros(Nperiod+1, NTrials);
i_curr = 0;
while i_curr < NTrials
    [x_smp,t_smp] = simByEuler(obj, Nperiod, 'DeltaTime', dt, 'NTrials', Nblock,...
        'Antithetic', Antithetic);
    p_smp(:, i_curr+ (1:Nblock)) = p(x_smp);
    i_curr = i_curr + Nblock;
end

save('twist_traj.mat', 'p_smp', 't_smp');
% x_smp = squeeze(x_smp);


epsilon_list = [0.5, 0.15, 0.1, 0.05];
max_smp = max(p_smp, [], 2);
q_smp = quantile(p_smp', 1-epsilon_list);

max_q = zeros(length(epsilon_list), 1);
for i = 1:length(max_q)
    max_q(i) = max(q_smp(i, :));
end
save('twist_traj_quantile.mat', 'max_q', 'epsilon_list', 'q_smp', 't_smp', 'max_smp')


%% plot
figure(3)
clf
hold on
plot(t_smp, q_smp)
plot(t_smp, max_smp, 'k')
xlim([0, 5])
title('Quantile Bounds for Fig 1', 'Fontsize', 16)
xlabel('t')
ylabel('VaR(x_3(t))')
legend({'\epsilon=0.5', '\epsilon=0.15', '\epsilon=0.1', '\epsilon=0.05', 'max'}, 'location', 'northeast')
