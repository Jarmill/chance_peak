
x = sdpvar(1,1);
y = sdpvar(1,1);
z = sdpvar(1,1);

% k = 0.5;
k = 4;


objective = z*k + x;

soc_vec = [1-y; 2*x; 2*z];
soc_bnd = 1 + y;

BOX = 2;
cons = [[x; y].^2 <= BOX^2; cone(soc_vec, soc_bnd)];

sol = optimize(cons, -objective, sdpsettings('solver', 'mosek'));

xr = value(x);
yr = value(y);
zr = value(z);
objr = value(objective);

sr = sqrt(yr-xr^2);

disp(zr - sr)