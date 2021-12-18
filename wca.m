function F = wca(dist,sigma,e)

% WCA potential given by the distance 'dist' to another object
% based on Lennard Jones. Used for the particles and the walls

if abs(dist) <= (2^(1/6))*sigma
    F = - 24*e*(((2*(sigma^12))/(dist^14)) - ((sigma^6)/(dist^8)))*dist;
else
    F = 0;
end

end