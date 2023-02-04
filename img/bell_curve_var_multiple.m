B = 4;
N = 201;

mu = 0;
sigma = 1;
epsilon = 0.05;


t = linspace(-B, B, N);

y = normpdf(t, mu, sigma);

teps = norminv(1-epsilon, mu, sigma);
tg = t(t>=teps);


%% multiple plots
mu_list = [mu-1.5, mu,  mu+0.7];
sigma_list = [1.75*sigma, sigma, sigma*0.4];
figure(2)
clf
tiledlayout(length(mu_list), 1)
c = linspecer(2);
ax = [];
for i = 1:length(mu_list)
    axcurr = nexttile;
    hold on
    y_curr = normpdf(t, mu_list(i), sigma_list(i));
    teps_curr = norminv(1-epsilon, mu_list(i), sigma_list(i));
%     plot([1,1]*teps_curr, ylim, '--k')
    
    if i == 2
        text(teps-0.9, 0.05, 'max VaR')
    else
        tg_curr = t(t>=teps_curr);
        patch([tg_curr, tg_curr(end:-1:1)], [y_curr(t>=teps_curr), zeros(1, length(tg_curr))],...
        c(2, :), 'edgecolor', 'none')
    end    
    patch([tg, tg(end:-1:1)], [y_curr(t>=teps), zeros(1, length(tg))],...
    'k', 'edgecolor', 'none')
plot(t, y_curr, 'LineWidth', 3, 'color', c(1, :));
    
    axis off
    ax = [ax; axcurr];   
end
linkaxes(ax, 'y')
