mpol('t');
mpol('x', 2);

supp = [t*(1-t)>=0; x'*x <= 4];

vars = struct('t', t, 'x', x);
% mb = meas_base(vars, supp);

% dynamics
f = [1 3; 4 -1]*x;
g = [1; x(2)+1];

dyn = struct('f', f, 'g', g);

d = 2;

% location support 
lsupp = chance_support(vars, 0.05);
% lsupp = lsupp.set_box(4);
lsupp = lsupp.set_box([-1, 3; -1.5, 2]);
% lsupp.X_init = X0;
lsupp.Tmax = 10;
% lsupp.mom_init = init_mom_handle;

Xsupp = lsupp.supp_sys();

SS_std = subsystem(lsupp, f);
SS = subsystem_sde(lsupp, dyn);

% mb = meas_base(vars, Xsupp);
% mb.mom_lie(d, [t; x], f)
% mb.mom_hess(d, [t; x], g)

Ay = SS.cons_liou(3)