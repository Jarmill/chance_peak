%Plot the pdf/likelihood of the gbm distribution
%Ito-closed-form solution given in https://en.wikipedia.org/wiki/Geometric_Brownian_motion
%parameters from example 1 of https://arxiv.org/pdf/2208.10752.pdf

b = -0.1;
sigma = sqrt(2)/2;

x0 = 0.25;
% T = 1;
% T = 5;
T = 10;
% Nt = 50;
Nt = 400;
t = linspace(0, T, Nt+1);
%the distribution at time 0 is a dirac-delta
t = t(2:end);

%fokker planck equation from https://www.mdpi.com/1099-4300/22/12/1432/htm
p = @(t, x) lognpdf(x/x0, b*t, sigma*sqrt(t));
Nx = 120;
Xmax = 5;

% x=linspace(0, Xmax, Nx);
x = logspace(-6, log10(Xmax), Nx);

pX = zeros(Nx, Nt);

for k = 1:Nt
    pX(:, k) = p(t(k), x);
end

[tt, xx] = meshgrid(t, x);

figure(3)
clf
subplot(1, 2, 1)
surf(tt, xx, log(pX));
zlim([-10, 1])

subplot(1, 2, 2)
surf(tt, xx, pX)
xlabel('t')
ylabel('x')
zlabel('p_t(x)')


% F = @(t,X) b * X;
% G = @(t,X) sigma * X;
% 
% x0 = 0.25;
% obj = sde(F, G, 'StartState', x0);    % dX = F(t,X)dt + G(t,X)dW
% 
% dt = 1e-3;
% T = 1;
% 
% %Options
% Nperiod = ceil(T/dt);
% Antithetic = true;
% NTrials = 1000;
% 
% [x_smp,t_smp] = simByEuler(obj, Nperiod, 'DeltaTime', dt, 'NTrials', NTrials,...
%     'Antithetic', Antithetic);
% x_smp = squeeze(x_smp);
% 
% %Plot
% figure(1)
% clf
% hold on
% plot(t_smp, x_smp);
% 
% titlestr = sprintf('dx = %0.2f x dt + %0.4f x dw', b, sigma);
% title(titlestr, 'fontsize', 16)
% xlabel('t')
% ylabel('x(t)')
% 
% xlim([0, T]);