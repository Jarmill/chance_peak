mset clear
clear all
%% variables

A = [-0.3, 0.8; -0.9, -0.1];

mpol('t', 1, 1)
mpol('x', 2, 1)
mpol('lam', 1, 1)

vars = struct;
vars.t = t;
vars.x = x;
vars.lam = lam;
Tmax = 10;


w0 = 0.25;
X0 = [-1; 0.5];


%% location support 

% lsupp = loc_support(vars);
lsupp = chance_support(vars);
lsupp = lsupp.set_box([-1, 1; -1, 1]*1.5);
lsupp.X_init = X0;
lsupp.lam_handle = @(d) NormalMom(d, d);
% lsupp.lam_handle = @(d) lam_max.^d;
lsupp.Tmax = 10;
lsupp.DISCRETE_TIME = 1;
lsupp.TIME_INDEP = 0;


%% testing peak estimation

%dynamics
% f = [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];

f = A*x - [0; 0.2*x(1)^2] + lam*w0*[prod(x); 0.1]; 

% f = x*(1-x);
% f = 0.8*x;
dyn = struct('f', f);

objective = -x(2);


% lsupp.epsilon =0.15;
lsupp.epsilon = 0.15;
lsupp.p_supp  = [-1.5; 1.5]; %max and min value of objective in X


% if SOLVE

order = 6;
PM = chance_peak_manager(lsupp, dyn, objective);
sol = PM.run(order, Tmax);
fprintf('VP at eps=%0.2f: %0.4f\n', lsupp.epsilon, sol.obj_rec)

epsilon_list = [0.15; 0.1; 0.05];
order_list = 1:6;
% peak_estimate = zeros(length(epsilon_list)+1, length(order_list));
% solver_time = zeros(length(epsilon_list)+1, length(order_list));

% start with the mean
for i = 1:length(order_list)
    lsupp.bound_type = 'mean';
    lsupp.epsilon = 0.5;
    PM = cvar_peak_manager(lsupp, dyn, objective);
    sol = PM.run(order_list(i), Tmax);
    peak_estimate(1, i) = sol.obj_rec;
    solver_time(1, i) = sol.solver_time;
    save('scatter_chance_test_sde_point_time.mat', 'epsilon_list', 'order_list', 'peak_estimate', 'solver_time');
end

% then do the chance-peak with a VP bound
for e = 1:length(epsilon_list)
% for e=1:1    
    lsupp.bound_type = 'vp';
    lsupp.epsilon = epsilon_list(e);
    for i = 1:length(order_list)
        PM = chance_peak_manager(lsupp, dyn, objective);
        sol = PM.run(order_list(i), Tmax);
        peak_estimate(e+1, i) = sol.obj_rec;
        solver_time(e+1, i) = sol.solver_time;
    end
    save('scatter_chance_test_sde_point_time.mat', 'epsilon_list', 'order_list', 'peak_estimate', 'solver_time');
end

