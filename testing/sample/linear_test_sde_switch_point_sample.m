mset clear
clear all
%% variables
SAMPLE = 1;
PLOT = 1;
PROCESS = 1;
mpol('t', 1, 1)
mpol('x', 2, 1)

vars = struct;
vars.t = t;
vars.x = x;

X0 = [0; 1];


%chance bound
epsilon = 0.1;
% epsilon = 0.05;


%% location support 

% lsupp = loc_support(vars);
lsupp = chance_support(vars, epsilon);
% lsupp = lsupp.set_box(4);
% lsupp = lsupp.set_box([-1, 3; -1.5, 2]);
lsupp = lsupp.set_box([-2, 2; -2, 2]);
lsupp.X_init = X0;
lsupp.Tmax = 5;


CHANCE = 1;
if CHANCE
%     lsupp.bound_type = 'cantelli';
    lsupp.bound_type = 'vp';
else
    lsupp.bound_type = 'mean';
end
%% testing peak estimation

%dynamics
f1 = [-5, -4; -1, -2]*x/2;
g1 = [0; 0.5*x(2)]/2;


f2 = [-2, -4; 5, -2]*x/2;
g2 = g1;
% f2 = [0.4 -0.5; 0.5 0.4]*x;
% g2 = [0.2; 0.2];

dyn = struct;
dyn.f = {f1, f2};
dyn.g = {g1, g2};

% dyn = struct('f', {f1, f2}, 'g', {g1, g2});

objective = x(2);

PM = chance_peak_manager(lsupp, dyn, objective);


%% Try out the SDE sampler
if SAMPLE
sampler = struct;
sampler.x = @() X0;

SMP = sampler_sde_uncertain(PM.loc, sampler);
SMP.mu = 0.3;

os_single = SMP.sample_traj(0, X0, [], 5);

Ntraj = 10;
osm = SMP.sample_traj_multi(Ntraj, lsupp.Tmax);
save('lin_test_switch_traj_big.mat', 'osm')
else
    load('lin_test_switch_traj_big.mat');
    Ntraj = length(osm);
end
%% plot the trajectory
if PLOT
figure(35)
N_highlight = 4;
c = linspecer(4);
clf
hold on
for i = 1:Ntraj

if i > (Ntraj - N_highlight)
    plot(osm{i}.x(:, 1), osm{i}.x(:, 2), 'color', c(Ntraj-i+1, :), 'LineWidth', 2);
else
    plot(osm{i}.x(:, 1), osm{i}.x(:, 2), 'c');
end
end
scatter(X0(1), X0(2), 200, 'ko')

plot(xlim, [1, 1]*[-0.335220126444672], 'r-.', 'LineWidth', 3)
plot(xlim, [1, 1]*[-0.654011542219790], 'k:', 'LineWidth', 3)
plot(xlim, [1, 1]*[-0.885312572563149], 'r', 'LineWidth', 3)


xlabel('$x_1$', 'interpreter', 'latex')
ylabel('$x_2$', 'interpreter', 'latex')
title('Stochastic Switched System', 'FontSize', 14)
end

%% get quantile bounds
if PROCESS
    % quantile bounds

epsilon_list = [0.5, 0.15, 0.1, 0.05, 0.01, 0.001]';
max_smp = max(p_smp, [], 2);
q_smp = quantile(p_smp', 1-epsilon_list);



% CVAR bounds

cvar_smp = zeros(size(q_smp));
for i = 1:size(cvar_smp, 1)
    for j = 1:size(cvar_smp, 2)
        pcurr = p_smp(j, :);
        cvar_smp(i, j) = mean(pcurr(pcurr>= q_smp(i, j)), 2);
    end    
end

% maximum outputs
max_q = zeros(length(epsilon_list), 1);
max_cvar = zeros(size(max_q));
for i = 1:length(max_q)
    max_q(i) = max(q_smp(i, :));
    max_cvar(i) = max(cvar_smp(i, :));
end

end
% pbaspect([diff(xlim), diff(ylim), 1])