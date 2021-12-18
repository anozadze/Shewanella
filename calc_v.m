function v = calc_v(vrange,npart)

% mean of the velocity distribution
mu_v = vrange/2;
% standard deviation of the velocity distribution
sigma_v = vrange/10;

% create a linspace array with npart+1 entries
v_dist = linspace(0,1,npart+2);
% crop 0 and 1 out of the array. This is necessary because of
% mathematical reasons (norminv of 0 and 1 is undefined!)
v_dist = v_dist([2:npart+1]);

% generate velocity distribution through the inversion method
v = norminv(v_dist,mu_v,sigma_v);

% randomly assigning the velocities to the particles
% while always taking the same seed, for reproducibility
rng('default');
rng(1);
v = v(randperm(length(v)));

end