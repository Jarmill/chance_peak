

mset clear
clear all
%% variables
SAMPLE = 0;
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

SOLVE = 1;

% if SOLVE
% PM = peak_manager(lsupp, f, objective);
PM = chance_peak_manager(lsupp, dyn, objective);


%% Try out the SDE sampler
if SAMPLE
sampler = struct;
sampler.x = @() X0;

SMP = sampler_sde_uncertain(PM.loc, sampler);
SMP.mu = 0.3;

os_single = SMP.sample_traj(0, X0, [], 5)

Ntraj = 300;
osm = SMP.sample_traj_multi(Ntraj, lsupp.Tmax);
save('lin_test_switch_traj.mat', 'osm')
else
    load('lin_test_switch_traj.mat');
    Ntraj = length(osm);
end
%% plot the trajectory
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

plot(xlim, [1, 1]*[-0.335220126444672], 'r--', 'LineWidth', 3)
plot(xlim, [1, 1]*[-0.736103602023845], 'r', 'LineWidth', 3)


xlabel('$x_1$', 'interpreter', 'latex')
ylabel('$x_2$', 'interpreter', 'latex')
title('Stochastic Switched System', 'FontSize', 14)

% pbaspect([diff(xlim), diff(ylim), 1])