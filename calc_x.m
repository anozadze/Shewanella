function x = calc_x(wall,npart)

% mean of the particle distribution
mu_x = wall/2;
% standard deviation of the particle distribution
sigma_x = wall/10;

% create a linspace array with npart+1 entries
x_dist = linspace(0,1,npart+2);
% crop 0 and 1 out of the array. This is necessary because of
% mathematical reasons (norminv of 0 and 1 is undefined!)
x_dist = x_dist([2:npart+1]);

% generate particle distribution through the inversion method
x = norminv(x_dist,mu_x,sigma_x);

end