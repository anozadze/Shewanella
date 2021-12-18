function F = force(pp,x,wall,sigma_pot,e)

% this function calculates the force of adjacent particles (for r_cutoff=1)
% and walls on particle pp.
%%note that the walls are treated like static particles here!

% reset the force (particles and walls)
F = 0.;

% the cutoff-radius determines which neighbours to look at
%%in 1D the cutoff radius is unnecessary.
%%In 2D it will become very important.
r_cutoff = 1;

% define a larger x that includes the walls (treated like particles!)
x_walls = [0 x wall];

% calculate the index of the next particle
%%in 1D the cutoff radius is unnecessary.
%%In 2D it will become very important.
pp_min = (pp+1) - r_cutoff;
pp_max = (pp+1) + r_cutoff;

% evaluate distance to neighboring particles (and walls) for force
%%with r_cutoff > 0, the vector is not possible anymore!
for fp = [pp_min pp_max]
    dx_part = x_walls(fp) - x(pp);
    F = F + wca(dx_part,sigma_pot,e);
end

end