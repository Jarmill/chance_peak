mset clear
clear all
SAVE = 0;
%% variables

mpol('t', 1, 1)
mpol('x', 2, 1)
mpol('y', 2, 1)

vars = struct;
vars.t = t;
vars.x = x;
vars.y = y;

X0 = [0; 0.75];
% X0 = [0; 1];

%% testing peak estimation

%dynamics
f1 = [-5, -4; -1, -2]*x/2;
% g1 = [0; 0.5*x(2)]/2;
g1 = [0; 0.1];


f2 = [-2, -4; 5, -2]*x/2;
g2 = g1;

dyn = struct;
% dyn.f = {f2};
% dyn.g = {g2};

dyn.f = f2;
dyn.g = g2;

%unsafe set
% theta_c = 5*pi/4; 
% Cu = [-1; -1];
% Cu = [ -0.5; -0.7];
% Ru = 0.5;
% Ru = 0.3;


% theta_c = 3*pi/2;
% theta_c = pi;

% 
% theta_c = 3*pi/2;
theta_c = 5*pi/4;
% Ru = 0.4;
Ru = 0.1;
% Cu = [-0.8; -0.8];
Cu = [-0.5;  -0.5];

% theta_c = 5*pi/4;
% Ru = 0.3;
% Cu = [-0.8; -0.8];
% Cu = [-0.8;  -0.5];
% theta_c = 
c1f = Ru^2 - (y(1) - Cu(1)).^2 - (y(2) - Cu(2)).^2;

w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(y(1) - Cu(1)) + w_c(2) * (y(2) - Cu(2)); 
% c2f = [];
unsafe_cons = [c1f; c2f];

% dyn = struct('f', {f1, f2}, 'g', {g1, g2});




%% location support 

lsupp = chance_distance_support(vars);
% lsupp = lsupp.set_box(4);
% lsupp = lsupp.set_box([-1, 3; -1.5, 2]);
box = [-0.75, 0.5; -1, 1];
lsupp = lsupp.set_box(box);
lsupp.X_init = X0;
lsupp.Tmax = 5;
lsupp.X_unsafe = unsafe_cons >= 0;
% lsupp.X_unsafe = (y==Cu);
lsupp.dist = (x-y)'*(x-y);
lsupp.epsilon = 0.15;

% objective = -x(1);

SOLVE = 1;

% if SOLVE
% PM = peak_manager(lsupp, f, objective);
% PM = chance_peak_manager(lsupp, dyn, objective);

% epsilon_list = [0.15; 0.1; 0.05];
epsilon_list = [0.15];
% order_list = 1:6;
% order_list = 4;
DO_MEAN = 1;
order_list = 1:5;
peak_estimate = zeros(length(epsilon_list)+1, length(order_list));
status = zeros(length(epsilon_list)+1, length(order_list));
solver_time = zeros(length(epsilon_list)+1, length(order_list));

if DO_MEAN
%start with the mean
for i = 1:length(order_list)
    lsupp.bound_type = 'mean';
    lsupp.epsilon = 0.5;
    PM = chance_distance_manager(lsupp, dyn);
    sol = PM.run(order_list(i));
    peak_estimate(1, i) = sol.obj_rec;
     status(1, i) = sol.status;
     solver_time(1, i) = sol.solver_time;
     if SAVE
    save('lin_test_dist_test_time_corr.mat', 'peak_estimate', 'status', 'order_list', 'epsilon_list', 'solver_time');
     end
end
end

disp(peak_estimate)

%then do the chance-peak with a VP bound
for e = 1:length(epsilon_list)
% for e=1:1    
    lsupp.bound_type = 'vp';
%      lsupp.bound_type = 'cantelli';
    lsupp.epsilon = epsilon_list(e);
    for i = 1:length(order_list)
        PM = chance_distance_manager(lsupp, dyn);
        sol = PM.run(order_list(i));
        peak_estimate(e+1, i) = sol.obj_rec;
         status(e+1, i) = sol.status;
         solver_time(e+1, i) = sol.solver_time;
         if SAVE
    save('lin_test_dist_test_time_corr.mat', 'peak_estimate','status',  'order_list', 'epsilon_list', 'solver_time');
         end
    end
end


