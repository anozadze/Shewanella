function [flag,dt] = flag_switch(pp,wall,x_check,x_new,dt)

% this function checks whether the particle has switched or crossed the wall
% which would mean that the particle pp is too fast!
% and then adjusts the time step and assigns the flag

% define a larger x that includes the walls (treated like particles!)
x_walls = [0 x_check wall];

% calculate the distance to the left and the right neighbour (or wall)
%%in 2D, left and right will not work anymore!
dx_left = x_new(pp) - x_walls((pp+1) - 1);
dx_right = x_walls((pp+1) + 1) - x_new(pp);

% check whether particle has switched with its neighbour or crossed the wall
if dx_left <= 0 || dx_right <= 0
    % decrease the time step, adjust the flag to 1
    dt = dt/2;
    flag = 1;
else
    % else leave the while loop by adjusting the flag to 0
    flag = 0;
end

end