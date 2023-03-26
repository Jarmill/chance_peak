%% data analysis

% load('motion_traj_small.mat')
load('motion_traj.mat')
p = @(x) -x(:, 2, :);

p_smp = squeeze(p(x_smp));

epsilon_list = [0.5, 0.15, 0.1, 0.05];
max_smp = max(p_smp, [], 2);
q_smp = quantile(p_smp', 1-epsilon_list);

max_q = zeros(length(epsilon_list), 1);
for i = 1:length(max_q)
    max_q(i) = max(q_smp(i, :));
end
save('flow_traj_quantile.mat', 'max_q', 'epsilon_list', 'q_smp', 't_smp', 'max_smp')

%% plot
figure(3)
clf
hold on
plot(t_smp, q_smp)
plot(t_smp, max_smp, 'k')
xlim([0, 5])
title('Quantile Bounds for Fig 1', 'Fontsize', 16)
xlabel('t')
ylabel('VaR(-x_2(t))')
legend({'\epsilon=0.5', '\epsilon=0.15', '\epsilon=0.1', '\epsilon=0.05', 'max'}, 'location', 'northeast')