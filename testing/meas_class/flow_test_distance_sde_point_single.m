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
X0 = [1; 1];

%chance bound
% epsilon = 0.1;
epsilon = 0.15;


%unsafe set
theta_c = 5*pi/4; 
Cu = [-0.5; -0.75];
Ru = 0.5;

c1f = Ru^2 - (y(1) - Cu(1)).^2 - (y(2) - Cu(2)).^2;

% theta_c = 3*pi/2;
% theta_c = 
w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(y(1) - Cu(1)) + w_c(2) * (y(2) - Cu(2)); 
unsafe_cons = [c1f; c2f];



%% location support 

% lsupp = loc_support(vars);
lsupp = chance_distance_support(vars, epsilon);
% lsupp = lsupp.set_box(4);
% lsupp = lsupp.set_box([-1, 3; -1.5, 2]);
lsupp = lsupp.set_box([-1.5, 1.5; -1.5, 1.5]);
lsupp.X_init = X0;
lsupp.Tmax = 5;
lsupp.X_unsafe = [c1f; c2f] >= 0;
lsupp.dist = (x-y)'*(x-y);
% lsupp.epsilon = ep;


CHANCE = 1;
if CHANCE
%     lsupp.bound_type = 'cantelli';
    lsupp.bound_type = 'vp';
else
    lsupp.bound_type = 'mean';
end
%% testing peak estimation

%dynamics
sigma = 0.1;
f = [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];
g = sigma*[0; 1];

dyn = struct('f', f, 'g', g);

% objective = -x(2);

SOLVE = 1;

PM = chance_distance_manager(lsupp, dyn);

%generate constraints
order = 4; %starting X0=C0, order 2: 0.5723, order 3: 0.5532
d = 2*order;
sol = PM.run(order, 5);
disp(sol.obj_rec)



