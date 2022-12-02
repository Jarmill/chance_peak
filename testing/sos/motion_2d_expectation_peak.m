%2d flow motion when drift term is constant
%find the maximal mean-expectation

% variables and dynamics
t = sdpvar(1,1);
x = sdpvar(2,1);

% b = -0.1;
sigma = 0.1;
f =  [x(2); -(x(1) +x(2) + 0.5*x(1)^3)];
g = sigma * [0;1];


p = -x(2); %objective


% order = 1; %-1.3504e-08
% order = 2; %0
% order = 3; %0
% order = 4; %0 
% order = 5; %0
order = 6; %0
d = 2*order;


x0 = [1; 1];

% k = sqrt(1/epsilon - 1); %cantelli bound
% k = sqrt(4/(9*epsilon) - 1); %VP bound

%% Support Sets
T = 5;
Xmax = 1.5;
Xall = struct('ineq', [t*(T-t); x.*(Xmax-x)], 'eq', []);

%% polynomials
%polynomial definition
[v, cv, mv] = polynomial([t;x], d);

Lv = jacobian(v, t) + jacobian(v, x)*f + 0.5*g'*hessian(v, x)*g;

v0 = replace(v, [t;x], [0; x0]);

objective = v0;

pcost = p;


[put_lie, conslie, coefflie] = constraint_psatz(-Lv, Xall, [t;x], d);
[put_cost, conscost, coeffcost] = constraint_psatz(v - pcost, Xall, [t;x], d);

cons = [conslie; conscost];
coeff = [coefflie; coeffcost; cv];

opts = sdpsettings('solver', 'mosek');


[sol,u,Q] = solvesos(cons,objective,opts,coeff);
% , X, vars, d)

%% recovery
phi_rec = value(phi);
lam_rec = value(lam);
v_rec = value(cv)'*mv;
v0_rec = replace(v_rec, [t; x], [0; x0]);
obj_rec = value(objective);

disp(obj_rec)