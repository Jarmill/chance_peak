mset clear
clear all
%% variables

mpol('t', 1, 1)
mpol('x', 2, 1)

vars = struct;
vars.t = t;
vars.x = x;

X0 = [1; 1];

%chance bound
% epsilon = 0.1;
epsilon = 0.05;


%% location support 

% lsupp = loc_support(vars);
lsupp = chance_support(vars, epsilon);
% lsupp = lsupp.set_box(4);
% lsupp = lsupp.set_box([-1, 3; -1.5, 2]);
lsupp = lsupp.set_box([-1.5, 1.5; -1.5, 1.5]);
lsupp.X_init = X0;
lsupp.Tmax = 5;


CHANCE = 1;
if CHANCE
%     lsupp.bound_type = 'cantelli';
    lsupp.bound_type = 'vp';
else
    lsupp.bound_type = 'mean';
end
%% testing peak estimation

%dynamics
f = [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];
g = [0; 0];

dyn = struct('f', f, 'g', g);

objective = -x(2);

SOLVE = 1;

% if SOLVE
% PM = peak_manager(lsupp, f, objective);
PM = chance_peak_manager(lsupp, dyn, objective);

%generate constraints
order = 5; %starting X0=C0, order 2: 0.5723, order 3: 0.5532
d = 2*order;
sol = PM.run(order, 10);
disp(sol.obj_rec)
% [obj_p, mom_con, supp_con] = PM.peak_cons(d);
% sol = PM.peak_solve(obj_p, mom_con,supp_con);
% end
% if REC
%     PM = PM.dual_process(d, sol.dual_rec);
% end


