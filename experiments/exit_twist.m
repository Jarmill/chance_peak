mset clear
clear all
%% variables

mpol('t', 1, 1)
mpol('x', 3, 1)

vars = struct;
vars.t = t;
vars.x = x;

X0 = [-0.5; 0; 0];

%% location support 

% lsupp = loc_support(vars);
lsupp = chance_support(vars);
lsupp.X = (sum(x.^4)<=1);
lsupp.X_term = (sum(x.^4)==1);
lsupp.X_init = X0;
lsupp.X_init = X0;
lsupp.Tmax = 5;

%% testing peak estimation

%dynamics
A_true = [-1 1 1; -1 0 -1; 0 1 -2];
B_true = [-1 0 -1;
          0 1 1;
          1 1 0]/2;

sigma = 0.1;
F = @(t,x) A_true*x - B_true*(4*x.^3 - 3*x);
G = @(t,x) sigma * [0; 0;1];

f = F(t, x);
g = G(t, x);


dyn = struct('f', f, 'g', g);


% objective = x(3);
objective = t;


% if SOLVE
% PM = peak_manager(lsupp, f, objective);
PM = chance_peak_manager(lsupp, dyn, objective);

epsilon_list = [0.15; 0.1; 0.05];
order_list = 1:5;
peak_estimate = zeros(length(epsilon_list)+1, length(order_list));


%start with the mean
for i = 1:length(order_list)
    lsupp.bound_type = 'mean';
    PM = chance_peak_manager(lsupp, dyn, objective);
    sol = PM.run(order_list(i));
    peak_estimate(1, i) = sol.obj_rec;
    save('exit_twist_test.mat', 'peak_estimate', 'order_list', 'epsilon_list');
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
        save('exit_twist_sde_test.mat', 'peak_estimate', 'order_list', 'epsilon_list');
    end
end

