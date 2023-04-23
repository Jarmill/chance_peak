%requires https://github.com/Jarmill/distance
mset clear
clear all
%% variables

mpol('t', 1, 1)
mpol('x', 2, 1)
mpol('y', 2, 1)

vars = struct;
vars.t = t;
vars.x = x;
vars.y = y;
X0 = [1; 0];
T = 5;


%unsafe set
theta_c = 5*pi/4; 
Cu = [-1; -1];
Ru = 0.5;

c1f = Ru^2 - (y(1) - Cu(1)).^2 - (y(2) - Cu(2)).^2;

% theta_c = 3*pi/2;
% theta_c = 
w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(y(1) - Cu(1)) + w_c(2) * (y(2) - Cu(2)); 
unsafe_cons = [c1f; c2f];

epsilon = 0.15;

%% location support 

% lsupp = loc_support(vars);
lsupp = chance_distance_support(vars, epsilon);
% lsupp = lsupp.set_box(4);
% lsupp = lsupp.set_box([-1, 3; -1.5, 2]);
box = [-1.25, 1.25; -1.5, 0.75];
lsupp = lsupp.set_box(box);
lsupp.X_init = X0;
lsupp.Tmax = 5;
lsupp.X_unsafe = unsafe_cons >= 0;
lsupp.dist = (x-y)'*(x-y);
lsupp.epsilon = 0.15;


CHANCE = 1;
if CHANCE
    lsupp.bound_type = 'cantelli';
%     lsupp.bound_type = 'vp';
else
    lsupp.bound_type = 'mean';
end
%% testing peak estimation

%dynamics
sigma = 0.1;
f = [x(2); -x(1) - (0.5).* x(1).^3 - x(2)];
g = sigma*[0; 1];

dyn = struct('f', f, 'g', g);

% objective = -x(2);

SOLVE = 1;

PM = chance_distance_manager(lsupp, dyn);

%% perform the experiment
if SOLVE
    
% epsilon_list = [0.15; 0.1; 0.05];
% order_list = 1:6;

epsilon_list = [0.15];
% order_list = 1:4;
% order_list = 5;
order_list = 7;
peak_estimate = zeros(length(epsilon_list)+1, length(order_list));
status = zeros(length(epsilon_list)+1, length(order_list));
solver_time = zeros(length(epsilon_list)+1, length(order_list));


%start with the mean
for i = 1:length(order_list)
    lsupp.bound_type = 'mean';
    lsupp.epsilon = 0.5;
    PM = chance_distance_manager(lsupp, dyn);
    sol = PM.run(order_list(i), T);
    peak_estimate(1, i) = sol.obj_rec;
    status(1, i) = sol.status;
    solver_time(1, i) = sol.solver_time;
    save('motion_distance_sde_test_time_7.mat', 'peak_estimate', 'order_list', 'epsilon_list', 'solver_time');
end

%then do the chance-peak with a VP bound
for e = 1:length(epsilon_list)
% for e=1:1    
%     lsupp.bound_type = 'cantelli';
    lsupp.bound_type = 'vp';
    lsupp.epsilon = epsilon_list(e);
    for i = 1:length(order_list)
       PM = chance_distance_manager(lsupp, dyn);
        sol = PM.run(order_list(i), T);
        peak_estimate(e+1, i) = sol.obj_rec;
        status(e+1, i) = sol.status;
        solver_time(e+1, i) = sol.solver_time;
        save('motion_distance_sde_test_time_7.mat', 'peak_estimate', 'order_list', 'epsilon_list', 'solver_time');
    end
end


end




