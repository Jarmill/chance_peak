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

box_setup = [-0.6, 0.6;
              -1, 1;
              -1, 1.5];
             
          
epsilon = 0.1;          
% lsupp = loc_support(vars);
lsupp = chance_support(vars, epsilon);
lsupp = lsupp.set_box(box_setup);
% lsupp = lsupp.set_box(4);
% lsupp = lsupp.set_box([-1, 3; -1.5, 2]);
% lsupp = lsupp.set_box([-1.5, 1.5; -1.5, 1.5]);
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

objective = x(3);

SOLVE = 1;

% if SOLVE
% PM = peak_manager(lsupp, f, objective);
PM = chance_peak_manager(lsupp, dyn, objective);

%generate constraints
order = 4; %starting X0=C0, order 2: 0.5723, order 3: 0.5532
d = 2*order;
sol = PM.run(order);
disp(sol.obj_rec)
% [obj_p, mom_con, supp_con] = PM.peak_cons(d);
% sol = PM.peak_solve(obj_p, mom_con,supp_con);
% end
% if REC
%     PM = PM.dual_process(d, sol.dual_rec);
% end


