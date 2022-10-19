%geometric brownian motion
%there is a semi-explicit expression for this
%the fokker-planck equation states that 

%% variables and dynamics
t = sdpvar(1,1);
x = sdpvar(1,1);

b = -0.1;
sigma = sqrt(2)/2;
f =  b * x;
g = sigma * x;

order = 2;
d = 2*order;


x0 = 0.25;

epsilon = 0.1;
p = x; %objective

%optimization variables
phi = sdpvar(3, 1);
lam = sdpvar(1, 1);


% k = sqrt(1/epsilon - 1); %cantelli bound
k = sqrt(4/(9*epsilon) - 1); %VP bound

%% Support Sets
T = 1;
Xmax = 2.5;
Xall = struct('ineq', [t*(T-t); x*(Xmax-x)], 'eq', []);

%% polynomials
%polynomial definition
[v, cv, mv] = polynomial([t,x], d);

Lv = jacobian(v, t) + jacobian(v, x)*f + g'*hessian(v, x)*g;

v0 = replace(v, [t;x], [0; x0]);

objective = v0 + lam - phi(1);

pcost = (1 - 2*phi(3))*p + (phi(1) + lam)*(p^2);


[plie, conslie, coefflie] = constraint_psatz(-Lv, Xall, [t;x], d);
[pcost, conscost, coeffcost] = constraint_psatz(v - pcost, Xall, [t;x], d);

cons = [cone(phi, lam); conslie; conscost];
coeff = [coefflie; coeffcost; phi; lam; cv];

opts = sdpsettings('solver', 'mosek');


[sol,u,Q] = solvesos(cons,objective,opts,coeff);
% , X, vars, d)
