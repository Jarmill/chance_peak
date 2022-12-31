mset clear
clear all
%% variables

mpol('t', 1, 1)
mpol('x', 2, 1)

vars = struct;
vars.t = t;
vars.x = x;


%do the box instead
% %initial set
C0 = [1; 1];
% C0 = 
% C0 = [0; 1.5];
% C0 = [1.5; ];
R0 = 0.4;

%% try and get Lebesgue moments of initial set


% y = LebesgueSphereMom(dv,r)
box = [1; 1]*C0' + [-1, -1; 1, 1]*R0/2;

smp = struct();
smp.x = @() (2*rand(2, 1)-1)*R0/2 + C0;

n=2;
init_mom_handle = @(d)  LebBoxMom_time_delta(d, box, 0, 0, 1);


%% location support 

lsupp = loc_support(vars);
% lsupp = chance_support(vars, 0.05);
% lsupp = lsupp.set_box(4);
% lsupp = lsupp.set_box([-3, 3; -3, 3]);
lsupp = lsupp.set_box([-1, 2; -1, 1.5]);
% lsupp = lsupp.set_box([-1, 3; -1.5, 2]);
% lsupp = lsupp.set_box([-1, 2; -1, 1.5]);
lsupp.Tmax = 10;
lsupp.mom_init = init_mom_handle;

%dynamics
f = [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];
objective = -x(2);

loc = location(lsupp, f, objective);


%% sampler
S = sampler_base(loc, smp);

N_sample = 150;
osm = S.sample_traj_multi(N_sample, lsupp.Tmax);

%% plot
figure(1)
clf
hold on
for i = 1:N_sample
    plot(osm{i}.x(:, 1), osm{i}.x(:, 2), 'c')
end
plot(box([1, 2, 2, 1, 1], 1), box([1, 1, 2, 2, 1], 2), 'k')

plot(xlim, [1, 1]*[-0.539607233605014], 'r--', 'LineWidth', 3)
plot(xlim, [1, 1]*-[0.936398391321781], 'r', 'LineWidth', 3)
xlabel('$x_1$', 'interpreter', 'latex')
ylabel('$x_2$', 'interpreter', 'latex')
title('Deterministic Flow System', 'FontSize', 14)
pbaspect([diff(xlim),diff(ylim),1])
