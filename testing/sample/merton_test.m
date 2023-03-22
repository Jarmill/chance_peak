rng(34, 'twister')
AssetPrice = 1;
Return = -0.1;
Sigma = 0.16;
% Sigma = sqrt(2)/2;
% JumpMean = 0.2;
JumpMean = 0.2;
JumpVol = 0.4;
JumpFreq = 4;
            
mertonObj = merton(Return,Sigma,JumpFreq,JumpMean,JumpVol,...
    'startstat',AssetPrice);
gbmObj = gbm(Return, Sigma, 'startstat', AssetPrice);

NTrials = 6;
dt = 1e-3;
T = 1;
Nperiods = ceil(T/dt);
[Paths, Times, Z] = mertonObj.simulate(Nperiods, 'DeltaTime', dt, 'NTrials', NTrials);
[PathG, TimeG] = gbmObj.simulate(Nperiods, 'DeltaTime', dt, 'NTrials', NTrials);
Paths = squeeze(Paths);
PathG = squeeze(PathG);


%% plot
figure(1)
clf
tiledlayout(1, 2)
a1 = nexttile;
plot(TimeG, PathG);
title('SDE (continuous)', 'fontsize', 16)
xlabel('t')
ylabel('x(t)')
box off
a2 = nexttile;
plot(Times, Paths);
linkaxes([a1; a2], 'xy')
xlim([0, T]);
xlabel('t')
ylabel('x(t)')
title('Levy Process (jumps)', 'fontsize', 16)
box off