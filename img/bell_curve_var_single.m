B = 4;
N = 201;

mu = 0;
sigma = 1;
epsilon = 0.05;


t = linspace(-B, B, N);

y = normpdf(t, mu, sigma);

teps = norminv(1-epsilon, mu, sigma);
tg = t(t>=teps);



%% single plot
figure(1)
clf
hold on
plot(t, y, 'LineWidth', 3);
patch([tg, tg(end:-1:1)], [y(t>=teps), zeros(1, length(tg))],...
    'k', 'edgecolor', 'none')
var_y = 0.04;
text(teps-0.5, var_y, 'VaR')
% ylim([var_y, max(y)])
axis off
