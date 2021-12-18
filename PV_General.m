%%for later: insert the documentation here
%

clear all

%% parameters for (repulsive) lennard-jones potential of particles, walls
sigma_pot = 0.5;
e = 10^-3;

%% parameters for execution

% var1: assign the values
niter = 2*10^3; %number of iterations %4.6*10^3
npart = 1*10^2; %number of particles
% implementing the simulation boundaries (in 1D: two walls, one at 0)
% for simplicity, walls are treated like static particles here.
wall = 5*10^2; %defines the simulation boundary
nplot = 1*10^2;

%% artefacts
% var2: ask for the values
% niter = input("\n number of iterations: ");
% npart = input("\n number of particles: ");
% implementing the simulation boundaries (in 1D: two walls, one at 0)
% for simplicity, walls are treated like static particles here.
% which means, they also exhibit a repulsive lennard-jones potential
% wall = 0;
%%

% particles should fit inside the walls! While also being properly distributed
%%the 8 is taken from trial and error experience. No logical reasons...
while wall < 8*npart*(2^(1/6))*sigma_pot
   wall = input("\n position of the simulation boundary (wall): ");
end

%% artefacts
%%alternatively:
%%walls = [0 rightwall];
%%nwalls = length(walls);
%%

%%artefacts
% to speed up the process
% nplot = input("\n plot after x iterations: ");
%%


% execution parameter part is finished.
% from here: better not change the code during testrun!!!

%%

% parameters for equation of motion
gamma = 1; % F - v*gamma = 0
dt0 = 0.01; % from experience: dt0 = 0.01 should work fine!
dt = dt0;

% defining the highest velocity of the particles
%%at the same time: speed of the fastest particle in the distribution
%%in reality: speed limit of the bacteria?
vrange = 20; %to enter manually? or to determine by trial and error?

% initializing new particle positions and directions for later
x_new = zeros(1,npart);
u_new = zeros(1,npart);

% calling functions for calculating particle position, direction, velocity
x = calc_x(wall,npart); % normal distribution
v = calc_v(vrange,npart); % normal distribution, but randomly assigned
u = calc_u(npart); % random, either 1 or -1

% for pdf plotting
fit_pd = fitdist(x','normal');
pd = normpdf(x,fit_pd.mu,fit_pd.sigma);
pd_max = max(pd);

% get the highest possible value for the velocity.
% Important for yrange when plotting
vmax = max(v);

%initializing particle positions for checking later
x_check = x;
%% artefacts
%v_check = v;
%u_check = u;
%%

% preparing for the first plot
% type helps to differentiate between initial and final state
type = 1;

% plotting the initial configuration of the particles
harryplotter(type,wall,pd_max,x,v,vmax,vrange,u,sigma_pot,npart)
pause(2) % pause to check initial configuration

% starting the main loop
for ii = 1:niter

    for pp = 1:npart

        % while loop for preventing switching of particles
        % flag for controlling the while loop

        flag = 1;
        while flag == 1

            x_new(pp) = 0;

            % calculate the forces on particle pp
            % due to the Lennard-Jones-Potential wca
            F = force(pp,x,wall,sigma_pot,e);

            % switch the direction in case of a repulsive force
            if (u(pp)*F < 0)
                u_new(pp) = -u(pp);
            else
                u_new(pp) = u(pp);
            end

            % calculate new position for particle pp
            dx = ((u_new(pp)*v(pp) + F/gamma))*dt;
            % dx = flag_switch_2(pp,npart,wall,x,x_new,dx);
            x_new(pp) = x(pp) + dx;

            % check whether the particle has switched with any other
            % or crossed one of the walls (treated like particles!)
            [flag,dt] = flag_switch(pp,wall,x_check,x_new,dt);

            %% artefacts
            % if it still doesn´t work...
            %v(pp) = v(pp)*dt/dt0;
            %%

        end

        % updating the positions for checking inside "flag_switch.m"
        % This is necessary since otherwise, the program doesn´t know the
        % which particles have already moved and which not. This leads to
        % errors in the whole program and might lead into an infinie loop!!
        x_check(pp) = x_new(pp);
        
        %% artefacts
        %v_check(pp) = v(pp);
        %u_check(pp) = u_new(pp);

        %x(pp) = x_new(pp);
        %u(pp) = u_new(pp);
        %%

        % reset the time step to the initial value
        % in case it was temporatily changed by "flag_switch.m"
        dt = dt0;

    end

    % updating the positions and directions
    x = x_new;
    u = u_new;
    
    % plot after every nplot steps
    if mod(ii,nplot) == 0
        % preparing for the real time plot
        % type helps to differentiate between initial and final state
        type = 2;

        % plotting the recent configuration of the particles
        harryplotter(type,wall,pd_max,x,v,vrange,vrange,u,sigma_pot,npart)
        %pause(2) % pause to check. Only implement if necessary!
    end

    %% artefacts
    %x = x_new;
    %v = v_new;
    %u = u_new;
    %%

end

