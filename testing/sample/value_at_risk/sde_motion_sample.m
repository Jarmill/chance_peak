rng(40, 'twister')
dt = 1e-3;
T =5;

% b = -0.1;
sigma = 0.1;
F = @(t,X) [X(2); -(X(1) +X(2) + 0.5*X(1)^3)];
% F = @(t,X) [X(2); -(X(1) +X(2) - (1/3)*X(1)^3)];
% G = @(t,X) sigma * [0; X(2)];
G = @(t,X) sigma * [0; 1];

p = @(x) squeeze(-x(:, 2, :));

x0 = [1; 1];
obj = sde(F, G, 'StartState', x0);    % dX = F(t,X)dt + G(t,X)dW

dt = 1e-3;
T =5;

Nperiod = ceil(T/dt);
Antithetic = true;
NTrials = 50000;
% NTrials = 4000;
Nblock = 1000;
p_smp = zeros(Nperiod+1, NTrials);
i_curr = 0;
while i_curr < NTrials
    [x_smp,t_smp] = simByEuler(obj, Nperiod, 'DeltaTime', dt, 'NTrials', Nblock,...
        'Antithetic', Antithetic);
    p_smp(:, i_curr+ (1:Nblock)) = p(x_smp);
    i_curr = i_curr + Nblock;
end

save('motion_traj.mat', 'p_smp', 't_smp');
% x_smp = squeeze(x_smp);

%% quantile bounds

epsilon_list = [0.5, 0.15, 0.1, 0.05, 0.01, 0.001];
max_smp = max(p_smp, [], 2);
q_smp = quantile(p_smp', 1-epsilon_list);



%% CVAR bounds

cvar_smp = zeros(size(q_smp));
for i = 1:size(cvar_smp, 1)
    for j = 1:size(cvar_smp, 2)
        pcurr = p_smp(j, :);
        cvar_smp(i, j) = mean(pcurr(pcurr>= q_smp(i, j)), 2);
    end    
end

%% maximum outputs
max_q = zeros(length(epsilon_list), 1);
max_cvar = zeros(size(max_q));
for i = 1:length(max_q)
    max_q(i) = max(q_smp(i, :));
    max_cvar(i) = max(cvar_smp(i, :));
end

%% save output
save('motion_traj_quantile_cvar.mat', 'max_cvar', 'max_q', 'epsilon_list', 'q_smp', 't_smp', 'max_smp')


% %% plot
% figure(3)
% clf
% hold on
% tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
% 
% 
% nexttile
% hold on
% plot(squeeze(x_smp1(:, 1, :)), squeeze(x_smp1(:, 2, :)), 'c')
% 
% xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 14)
% ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 14)
% title('Sample Paths', 'Fontsize', 16)
% % pbaspect([diff(xlim), diff(ylim),1])
% 
% cl = linspecer(4);
% for i = 1:4
%     plot(squeeze(x_smp1(:, 1, i)), squeeze(x_smp1(:, 2, i)), 'color', cl(i, :), 'linewidth', 3)
% end
% scatter(x0(1), x0(2), 200, 'k')
% % title('Value-at-Risk Bounds', 'Fontsize', 16)
% 
% nexttile;
% plot(t_smp, q_smp, 'linewidth', 3)
% % plot(t_smp, max_smp, 'k')
% xlim([0, 5])
% title('Value-at-Risk Bounds', 'Fontsize', 16)
% xlabel('$t$', 'interpreter', 'latex', 'fontsize', 14)
% ylabel('$VaR(p_\# \mu_{t})$', 'interpreter', 'latex', 'fontsize', 14)
% % legend({'\epsilon=0.5', '\epsilon=0.15', '\epsilon=0.1', '\epsilon=0.05', 'max'}, 'location', 'northeast')
% legend({'\epsilon=0.5', '\epsilon=0.15', '\epsilon=0.1', '\epsilon=0.05', '\epsilon=0.01', '\epsilon=0.001',}, ...
%     'location', 'northeast', 'fontsize', 12)
% 
% % figure(4)
% % clf
