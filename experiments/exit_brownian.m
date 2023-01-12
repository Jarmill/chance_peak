mset clear
clear all
%% variables

mpol('t', 1, 1)
mpol('x', 3, 1)

vars = struct;
vars.t = t;
vars.x = x;

X0 = [0; 0; 0];

%% location support 
% lsupp = loc_support(vars);
lsupp = chance_support(vars);
% lsupp = lsupp.set_box(4);
% lsupp = lsupp.set_box([-1, 3; -1.5, 2]);
% lsupp = lsupp.set_box(box_setup);
lsupp.X = (sum(x.^4)<=1);
lsupp.X_term = (sum(x.^4)==1);
% lsupp.X = (sum(x.^2)<=1);
% lsupp.X_term = (sum(x.^2)==1);
lsupp.X_init = X0;
lsupp.Tmax = 3;

%% testing peak estimation

%dynamics
A_true = [-1 1 1; -1 0 -1; 0 1 -2];
B_true = [-1 0 -1;
          0 1 1;
          1 1 0]/2;

sigma =1;
F = @(t,x) zeros(3, size(x, 2));
% F = @(t, x)  [1.1, -1, 0; 1, -1, -1; -1 1 -1]*x;
% G = @(t,x) sigma*ones(3, size(x, 2));
G = @(t, x) sigma*ones(3, size(x, 2));

f = F(t, x);
g = G(t, x);


dyn = struct('f', f, 'g', g);


% objective = sum(x.^2);
objective = t;


% if SOLVE
% PM = peak_manager(lsupp, f, objective);
PM = chance_peak_manager(lsupp, dyn, objective);

epsilon_list = [0.15; 0.1; 0.05];
order_list = 1:5;
peak_estimate = zeros(length(epsilon_list)+1, length(order_list));
status = zeros(length(epsilon_list)+1, length(order_list));

%start with the mean
for i = 1:length(order_list)
    lsupp.bound_type = 'mean';
    PM = chance_peak_manager(lsupp, dyn, objective);
    sol = PM.run(order_list(i));
    peak_estimate(1, i) = sol.obj_rec;
    status(1, i) = sol.status;
    save('exit_brownian_sde_test.mat', 'peak_estimate', 'order_list', 'epsilon_list');
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
        status(e+1, i) = sol.status;
        save('exit_brownian_sde_test.mat', 'peak_estimate', 'order_list', 'epsilon_list');
    end
end

