rng(40, 'twister')
dt = 1e-3;
T =5;


sigma = 0.1;


F = @(t,x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];
% F = @(t,X) [X(2); -(X(1) +X(2) - (1/3)*X(1)^3)];
% G = @(t,X) sigma * [0; X(2)];
% G = @(t,x) sigma * [0; 0; 1];
G = @(t,x) sigma * [0; 1];


box = [-1.1, 1.75; -1.5, 1.5];

be = @(t, x) box_event_nan(t, x, box);
% pe = @(t, e) = x*be(t, x) + NaN*ones(size(x))*(1-be(t, x));

% x0 = [1; 1];
x0 = [1; 1];
obj = sde(F, G, 'StartState', x0);    % dX = F(t,X)dt + G(t,X)dW

%unsafe set 
theta_c = 5*pi/4; 
Cu = [-0.5; -0.75];
Ru = 0.5;

p = @(x) -aff_half_circ_dist(x, Ru, theta_c, Cu)^2;

Nperiod = ceil(T/dt);
Antithetic = true;
NTrials = 50000;
Nblock = 1000;

% NTrials = 200;
% Nblock = 100;
p_smp = zeros(Nperiod+1, NTrials);
i_curr = 0;
while i_curr < NTrials
    [x_smp,t_smp] = simByEuler(obj, Nperiod, 'DeltaTime', dt, 'NTrials', Nblock,...
        'Antithetic', Antithetic, 'processes', be);   
    for i = 1:size(x_smp, 1)
        for j = 1:size(x_smp, 3)
            p_smp(i, i_curr+j) = p(squeeze(x_smp(i, :, j)));
        end
    end
%     p_smp(:, i_curr+ (1:Nblock)) = p(x_smp);
    i_curr = i_curr + Nblock;
end

save('flow_traj.mat', 'p_smp', 't_smp');
% x_smp = squeeze(x_smp);


epsilon_list = [0.5, 0.15, 0.1, 0.05];
max_smp = max(p_smp, [], 2);
q_smp = quantile(p_smp', 1-epsilon_list);

max_q = zeros(length(epsilon_list), 1);
for i = 1:length(max_q)
    max_q(i) = max(q_smp(i, :));
end
save('flow_traj_quantile.mat', 'max_q', 'epsilon_list', 'q_smp', 't_smp', 'max_smp')


%% plot
figure(3)
clf
hold on
plot(t_smp, q_smp)
plot(t_smp, max_smp, 'k')
xlim([0, 5])
title('Quantile Bounds for Fig 1', 'Fontsize', 16)
xlabel('t')
ylabel('VaR(-x_2(t))')
legend({'\epsilon=0.5', '\epsilon=0.15', '\epsilon=0.1', '\epsilon=0.05', 'max'}, 'location', 'northeast')


figure(4)
clf
plot(squeeze(x_smp(:, 1, :)), squeeze(x_smp(:, 2, :)), 'c')

theta_half_range = linspace(theta_c-pi/2, theta_c + pi/2, 200);
circ_half = [cos(theta_half_range); sin(theta_half_range)];
Xu = Cu + circ_half* Ru;
patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none', 'DisplayName', 'Unsafe Set')
    
  

%% functions

function ba = box_event_nan(t, x, box)
%     [be, term, dir] = box_event(t, x, box);
    
    be = all((x>= box(:, 1)) & (x <= box(:, 2)));

    if be
        ba = x;
    else
        ba = NaN(size(x));
    end

end

function dist_out = half_circ_dist(x_in, R)
    %return the L2 distance  between the point x_in and the half circle
    %||x_in||^2 <= R^2 intersect x_in(2) <= 0.
%     reshape(x_in, [], 1);
    if x_in(2) >= 0
        %flat region
        if x_in(1) < -R
            dist_out = hypot(x_in(1)+R, x_in(2));
        elseif x_in(1) > R
            dist_out = hypot(x_in(1)-R, x_in(2));
        else
            dist_out = x_in(2);
        end
    else
        %circle region
        dist_out = max(norm(x_in, 2)-R, 0);
    end

end

function dist_out = aff_half_circ_dist(x_in, R, theta_c, Cu)

    theta_cf = theta_c - 3*pi/2;
    Rot_mat = [cos(theta_cf) -sin(theta_cf); sin(theta_cf) cos(theta_cf)];
    x_aff = Rot_mat'*(x_in - Cu);
    
    dist_out = half_circ_dist(x_aff, R);

end


function x_dist = dist_contour(Ntheta, R, c)
    %compute a contour at distance c away from the half-circle with N_theta
    %sample points


    theta_q1 = linspace(0, pi/2, Ntheta);
    theta_q2 = linspace(pi/2, pi, Ntheta);
    theta_q34 = linspace(pi, 2*pi, 2*Ntheta);

    %contour level
    

    x_right = [c*cos(theta_q1)+R; c*sin(theta_q1)];
    x_left = [c*cos(theta_q2)-R; c*sin(theta_q2)];
    x_bottom = [(c+R)*cos(theta_q34); (c+R)*sin(theta_q34)];
    x_dist = [x_right, x_left, x_bottom];
end
