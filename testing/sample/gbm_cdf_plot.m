b = -0.1;
sigma = sqrt(2)/2;

x0 = 0.25;
T = 1;
% T = 5;

% t = linspace(0, T, Nt+1);
%the distribution at time 0 is a dirac-delta
% t = t(2:end);

Nx = 200;
Xmax = 3;

% x=linspace(0, Xmax, Nx);
x = logspace(-5, log10(Xmax), Nx);

%fokker planck equation from https://www.mdpi.com/1099-4300/22/12/1432/htm
pcdf = @(x) logncdf(x/x0, b*T, sigma*sqrt(T));

figure(10)
semilogx(x, pcdf(x))