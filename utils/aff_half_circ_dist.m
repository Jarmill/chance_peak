function dist_out = aff_half_circ_dist(x_in, R, theta_c, Cu)


%find the distance between a point x_in and a half-circle with center Cu,
%radius R, and angle theta_c

Npt = size(x_in, 2);
dist_out = zeros(1, Npt);

theta_cf = theta_c - 3*pi/2;
Rot_mat = [cos(theta_cf) -sin(theta_cf); sin(theta_cf) cos(theta_cf)];

%iterate through the points
for i = 1:Npt    
    x_aff = Rot_mat'*(x_in(:, i) - Cu);
    
    dist_out(i) = half_circ_dist(x_aff, R);
end

end
