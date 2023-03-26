mset clear
clear all
%% variables

mpol('t', 1, 1)
mpol('x', 2, 1)

vars = struct;
vars.t = t;
vars.x = x;

X0 = [1; 1];


%% location support 

% lsupp = loc_support(vars);
lsupp = chance_support(vars);
% lsupp = lsupp.set_box(4);
% lsupp = lsupp.set_box([-1, 3; -1.5, 2]);
% lsupp = lsupp.set_box([-1, 1.5; -1.5, 1.5]);
lsupp = lsupp.set_box([-1, 1.5; -1.5, 1.5]);
lsupp.X_init = X0;
lsupp.Tmax = 5;

%% testing peak estimation

%dynamics
% f = [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];
f = [x(2); -x(1) - (0.5).* x(1).^3 - x(2)];
sigma = 0.1;
g = [0; sigma];

dyn = struct('f', f, 'g', g);

objective = -x(2);


% if SOLVE
% PM = peak_manager(lsupp, f, objective);
PM = chance_peak_manager(lsupp, dyn, objective);

epsilon_list = [0.15; 0.1; 0.05];
order_list = 1:6;

% epsilon_list = [0.15; 0.1; 0.05];
% order_list = 4;
peak_estimate = zeros(length(epsilon_list)+1, length(order_list));
solver_time = zeros(length(epsilon_list)+1, length(order_list));

%start with the mean
for i = 1:length(order_list)
    lsupp.bound_type = 'mean';
    PM = chance_peak_manager(lsupp, dyn, objective);
    sol = PM.run(order_list(i));
    peak_estimate(1, i) = sol.obj_rec;
    solver_time(1, i) = sol.solver_time;
    save('motion_test_sde_point_time.mat', 'epsilon_list', 'order_list', 'peak_estimate', 'solver_time');
end

%then do the chance-peak with a VP bound
for e = 1:length(epsilon_list)
% for e=1:1    
    lsupp.bound_type = 'vp';
    lsupp.epsilon = epsilon_list(e);
    for i = 1:length(order_list)
        PM = chance_peak_manager(lsupp, dyn, objective);
        sol = PM.run(order_list(i));
        peak_estimate(e+1, i) = sol.obj_rec;
        solver_time(e+1, i) = sol.solver_time;
    end
    save('motion_test_sde_point_time.mat', 'epsilon_list', 'order_list', 'peak_estimate', 'solver_time');
end

