mset clear
clear all
%% variables

mpol('t', 1, 1)
mpol('x', 2, 1)

vars = struct;
vars.t = t;
vars.x = x;

%% try and get Lebesgue moments of initial set
%form a box
C0 = [1; 1];
% C0 = [0; 1.5];
R0 = 0.4;
box = [1; 1]*C0' + [-1, -1; 1, 1]*R0/2;

n=2;
% dv = genPowGlopti(n,d);
% init_mom = LebesgueBoxMom( d, box, 0 );
init_mom_handle = @(d)  LebBoxMom_time_delta(d, box, 0, 0, 1);


%% location support 

% lsupp = loc_support(vars);
lsupp = chance_support(vars);
% lsupp = lsupp.set_box(4);
% lsupp = lsupp.set_box([-1.5, 1.5; -1.5, 1.5]);
lsupp = lsupp.set_box([-1, 2; -1, 1.5]);
% lsupp.X_init = X0;
lsupp.X_init = (x-box(1, :)').*(box(2, :)'-x) >=0 ;
lsupp.Tmax = 5;
% lsupp.mom_init = init_mom_handle;

%% parameters of dynamics
%dynamics
f = [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];
g = [0; 0];

dyn = struct('f', f, 'g', g);

objective = -x(2);

%% testing peak estimation


order_list = 1:6;
peak_estimate = zeros(1, length(order_list));




%start with the mean
for i = 1:length(order_list)
    lsupp.bound_type = 'mean';
    PM = chance_peak_manager(lsupp, dyn, objective);
    sol = PM.run(order_list(i));
    peak_estimate(1, i) = sol.obj_rec;
end


%% sample


% if SOLVE
% PM = chance_peak_manager(lsupp, dyn, objective);

% %generate constraints
% order = 3; 
% d = 2*order;
% sol = PM.run(order);
% disp(sol.obj_rec)


