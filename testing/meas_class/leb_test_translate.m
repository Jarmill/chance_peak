mset clear
clear all
%% variables

mpol('t', 1, 1)
mpol('x', 2, 1)
mpol('y', 2, 1)
d = 4;

vars = struct;
vars.t = t;
vars.x = x;

%initial set
C0 = [1.5; 0];
R0 = 0.4;


xc = (R0*x + C0);

dv = genPowGlopti(2,d);
yc = prod([xc'].^dv, 2); %initial measure
yunit = LebesgueSphereMom(dv,1);