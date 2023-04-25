mset clear
clear all
SAVE = 0;
%% variables

mpol('t', 1, 1)
mpol('x', 2, 1)

vars = struct;
vars.t = t;
vars.x = x;

X0 = [0; 0.75];

%% location support 

% lsupp = loc_support(vars);
lsupp = chance_support(vars);
box = [-1.25, 1; -1, 1];
lsupp = lsupp.set_box(box);
% lsupp = lsupp.set_box([-2, 2; -2, 2]);
lsupp.X_init = X0;
lsupp.Tmax = 5;

%% testing peak estimation

%dynamics
f1 = [-5, -4; -1, -2]*x/2;
g1 = [0; 0.1];


f2 = [-2, -4; 5, -2]*x/2;
g2 = g1;

dyn = struct;
dyn.f = {f2};
dyn.g = {g2};

% dyn = struct('f', {f1, f2}, 'g', {g1, g2});

% objective = -x(2);

Cu = [-1;  -1];
objective = -sum((x-Cu).^2);
% objective = -x(1);

SOLVE = 1;

% if SOLVE
% PM = peak_manager(lsupp, f, objective);
PM = chance_peak_manager(lsupp, dyn, objective);

% epsilon_list = [0.15; 0.1; 0.05];
% epsilon_list = [0.15; 0.1; 0.05];
epsilon_list = [0.15];
% order_list = 1:6;
order_list = 1:5;
% order_list = 1:3;
peak_estimate = zeros(length(epsilon_list)+1, length(order_list));
status = zeros(length(epsilon_list)+1, length(order_list));
solver_time = zeros(length(epsilon_list)+1, length(order_list));


%start with the mean
for i = 1:length(order_list)
    lsupp.bound_type = 'mean';
    lsupp.epsilon = 0.5;
    PM = chance_peak_manager(lsupp, dyn, objective);
    sol = PM.run(order_list(i));
    peak_estimate(1, i) = sol.obj_rec;
     status(1, i) = sol.status;
     solver_time(1, i) = sol.solver_time;
     if SAVE
    save('lin_test_circle_test_time_corr.mat', 'peak_estimate', 'status', 'order_list', 'epsilon_list', 'solver_time');
     end
end

disp(peak_estimate)

%then do the chance-peak with a VP bound
for e = 1:length(epsilon_list)
% for e=1:1    
    lsupp.bound_type = 'vp';
%     lsupp.bound_type = 'cantelli';
    lsupp.epsilon = epsilon_list(e);
    for i = 1:length(order_list)
        PM = chance_peak_manager(lsupp, dyn, objective);
        sol = PM.run(order_list(i));
        peak_estimate(e+1, i) = sol.obj_rec;
         status(e+1, i) = sol.status;
         solver_time(e+1, i) = sol.solver_time;
         if SAVE
    save('lin_test_circle_test_time_corr.mat', 'peak_estimate','status',  'order_list', 'epsilon_list', 'solver_time');
         end
    end
end


