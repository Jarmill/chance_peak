mset clear
clear all
%% variables

mpol('t', 1, 1)
mpol('x', 2, 1)

vars = struct;
vars.t = t;
vars.x = x;


%do the box instead
% %initial set
C0 = [1.5; 0];
R0 = 0.4;

% X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);

X0 = C0;
% X0 = [];

%% try and get Lebesgue moments of initial set


% y = LebesgueSphereMom(dv,r)
box = [1; 1]*C0' + [-1, -1; 1, 1]*R0/2;

n=2;
% dv = genPowGlopti(n,d);
% init_mom = LebesgueBoxMom( d, box, 0 );
init_mom_handle = @(d)  LebBoxMom_time_delta(d, box, 0, 0, 1);




%% location support 

% lsupp = loc_support(vars);
lsupp = chance_support(vars, 0.05);
% lsupp = lsupp.set_box(4);
lsupp = lsupp.set_box([-1, 3; -1.5, 2]);
% lsupp.X_init = X0;
lsupp.Tmax = 10;
lsupp.mom_init = init_mom_handle;
lsupp.bound_type = 'mean';
%% testing peak estimation

%dynamics
f = [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];
g = [0; 0];

dyn = struct('f', f, 'g', g);

objective = -x(2);

SOLVE = 1;
REC = 1;
% if SOLVE
% PM = peak_manager(lsupp, f, objective);
PM = chance_peak_manager(lsupp, dyn, objective);

%generate constraints
order = 2; %starting X0=C0, order 2: 0.5723, order 3: 0.5532
d = 2*order;
sol = PM.run(order, 10);
disp(sol.obj_rec)
% [obj_p, mom_con, supp_con] = PM.peak_cons(d);
% sol = PM.peak_solve(obj_p, mom_con,supp_con);
% end
% if REC
%     PM = PM.dual_process(d, sol.dual_rec);
% end


