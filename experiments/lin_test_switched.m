mset clear
clear all
SAVE = 0;
%% variables

mpol('t', 1, 1)
mpol('x', 2, 1)

vars = struct;
vars.t = t;
vars.x = x;

X0 = [0; 1];

%% location support 

% lsupp = loc_support(vars);
lsupp = chance_support(vars);
lsupp = lsupp.set_box([-2, 2; -2, 2]);
lsupp.X_init = X0;
lsupp.Tmax = 5;

%% testing peak estimation

%dynamics
f1 = [-5, -4; -1, -2]*x/2;
g1 = [0; 0.5*x(2)]/2;


f2 = [-2, -4; 5, -2]*x/2;
g2 = g1;

dyn = struct;
dyn.f = {f1, f2};
dyn.g = {g1, g2};

% dyn = struct('f', {f1, f2}, 'g', {g1, g2});

objective = -x(2);

SOLVE = 1;

% if SOLVE
% PM = peak_manager(lsupp, f, objective);
PM = chance_peak_manager(lsupp, dyn, objective);

epsilon_list = [0.15; 0.1; 0.05];
% epsilon_list = [0.15];
% order_list = 1:6;
order_list = 1:3;
peak_estimate = zeros(length(epsilon_list)+1, length(order_list));
status = zeros(length(epsilon_list)+1, length(order_list));


%start with the mean
for i = 1:length(order_list)
    lsupp.bound_type = 'mean';
    PM = chance_peak_manager(lsupp, dyn, objective);
    sol = PM.run(order_list(i));
    peak_estimate(1, i) = sol.obj_rec;
     status(1, i) = sol.status;
     if SAVE
    save('lin_test_switched_test.mat', 'peak_estimate', 'status', 'order_list', 'epsilon_list');
     end
end

disp(peak_estimate)

% %then do the chance-peak with a VP bound
% for e = 1:length(epsilon_list)
% % for e=1:1    
%     lsupp.bound_type = 'vp';
% %     lsupp.bound_type = 'cantelli';
%     lsupp.epsilon = epsilon_list(e);
%     for i = 1:length(order_list)
%         PM = chance_peak_manager(lsupp, dyn, objective);
%         sol = PM.run(order_list(i));
%         peak_estimate(e+1, i) = sol.obj_rec;
%          status(e+1, i) = sol.status;
%          if SAVE
%     save('lin_test_switched_test.mat', 'peak_estimate','status',  'order_list', 'epsilon_list');
%          end
%     end
% end


