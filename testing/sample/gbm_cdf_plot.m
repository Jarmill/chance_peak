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

%% CDF in Space
% x=linspace(0, Xmax, Nx);
x = logspace(-5, log10(Xmax), Nx);

%fokker planck equation from https://www.mdpi.com/1099-4300/22/12/1432/htm
pcdf = @(x) logncdf(x/x0, b*T, sigma*sqrt(T));

figure(10)
semilogx(x, pcdf(x))
xlabel('x')
ylabel('cdf(x|T)')
title('CDF evaluation at terminal time')

%% CDF in Time
figure(11);
%compare the extracted bound from geometric_1d_sos.m vs the true cdf in teh
%time interval
epsilon = 0.1;

x_5_1 = 0.765525507142643;  %order 5, Xmax = 3, T = 1
x_5_5 = 1.088064660165312; %order 5, Xmax = 3, T = 5
x_5_all_5 = 1.111474751861115; %order 5, Xmax = inf, T=5

x_5_all_Tscale = 1.300065273747971;
% x_5_cant = 1.706856693288101; %cantelli bound
pcdf_5 = @(t) logncdf(x_5_1/x0, b*t, sigma*sqrt(t)); 
t = linspace(0, T, Nx+1);
t = t(2:end);

% pt = zeros(size(t));
pt = pcdf_5(t);
% for i = 1:Nx
% %     pt(i)
%     pt_curr = pcdf_5(t(i));
%     pt(i) = pt_curr;
% end

clf
hold on
plot(t, pt)
plot([0, T], (1-epsilon)*[1,1], 'k:', 'LineWidth', 2);
hold off
xlabel('t')
ylabel('cdf(t | x_5)')
title(sprintf('CDF evaluation at order-5 certified epsilon=%0.2f bound', epsilon)) 
ylim([0, 1])