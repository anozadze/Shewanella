function c_local = local_c(xp,c0,x0)

% very simple function. Just gives back the local concentration c_local at
% the particle position xp (in the main program: "x(pp)") of particle pp.
% however, this might become a bit more complicated in 2D.

%% artefacts
%slowing down
%input: no vectors!
%avoid conflicting
%%

% exponential gradient with maximum at (0, c0)
c_local = c0*exp(-xp/x0);

end