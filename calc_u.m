function u = calc_u(npart)
    
    % initializing orientation (unit vector. In 1D: +1 or -1)
    u = zeros(1,npart);
    % assigning a random orientation (either +1 or -1) for each particle
    % while always taking the same seed, for reproducibility
    rng('default');
    rng(1);
    r_numbers = rand(1,npart);
    for iu = 1:length(u)
        rng('default');
        rng(1);
        u(iu) = 2*round(r_numbers(iu))-1;
        %%alternative: take cdf!!
    end
    
end