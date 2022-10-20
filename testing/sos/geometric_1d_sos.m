%geometric brownian motion
%there is a semi-explicit expression for this
%the fokker-planck equation states that 

% variables and dynamics
t = sdpvar(1,1);
x = sdpvar(1,1);

b = -0.1;
sigma = sqrt(2)/2;
f =  b * x;
g = sigma * x;


x0 = 0.25;

epsilon = 0.1;
p = x; %objective

%T = 1, Xmax = 3;
% order = 1;  %  1.6897
% order = 2; %   0.7850
% order = 3; % 0.7737
% order = 4; %0.7674
order = 5; %0.7655
% order = 6; %0.7621 (unknown)
% order = 7; %  0.7537 (unknown)
% order = 8; %    0.7469 (unknown)




d = 2*order;

%optimization variables
phi = sdpvar(3, 1);
lam = sdpvar(1, 1);


k = sqrt(1/epsilon - 1); %cantelli bound
% k = sqrt(4/(9*epsilon) - 1); %VP bound

%% Support Sets
% T = 1;
T = 5;
% Xmax = 3;
Xmax = 5;
% Xall = struct('ineq', [t*(T-t); x*(Xmax-x)], 'eq', []);
Xall = struct('ineq', [t*(T-t); x], 'eq', []);

%% polynomials
%polynomial definition
[v, cv, mv] = polynomial([t,x], d);

Lv = jacobian(v, t) + jacobian(v, x)*f + g'*hessian(v, x)*g;

v0 = replace(v, [t;x], [0; x0]);

objective = v0 + lam - phi(1);

pcost = (1 - 2*phi(3))*p + (phi(1) + lam)*(p^2);


[plie, conslie, coefflie] = constraint_psatz(-Lv, Xall, [t;x], d);
[pcost, conscost, coeffcost] = constraint_psatz(v - pcost, Xall, [t;x], d);

cons = [cone(phi, lam); conslie; conscost; phi(2)==(k/2)];
coeff = [coefflie; coeffcost; phi; lam; cv];

opts = sdpsettings('solver', 'mosek');


[sol,u,Q] = solvesos(cons,objective,opts,coeff);
% , X, vars, d)

%% recovery
phi_rec = value(phi);
lam_rec = value(lam);
v_rec = value(cv)'*mv;
v0_rec = replace(v_rec, [t; x], [0, x0]);
obj_rec = value(objective)