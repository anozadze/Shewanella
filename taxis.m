function up_new = taxis(c_local,c0,dc,up_new,seed)

% default turning probability
% value of 0.001 taken from master thesis
p0 = 0.001;

% adjustment parameters for the turning probability
% values of 10 taken from master thesis
k1 = 10; % right side
k2 = 10; % left side

c_idealmin = 0.4*c0; % left edge of the ideal concentration region
c_idealmax = 0.6*c0; % right edge of the ideal concentration region
c_buffermin = 0.3*c0; % left buffer
c_buffermax = 0.8*c0; % right buffer

% now checking where the particle is and assigning the right turning
% probability p_switch
if c_local > c_idealmin && c_local < c_idealmax
    % particle is inside the both edges (ideal concentration)
    p_turn = p0;
    %
elseif c_local > c_idealmax && c_local < c_buffermax && up_new*dc <= 0
    % particle is between right edge and right buffer
    % but right direction
    p_turn = p0;
    %
elseif c_local > c_buffermin && c_local < c_idealmin && up_new*dc >= 0
    % particle is between right edge and right buffer
    % but right direction
    p_turn = p0;
    %
elseif c_local > c_buffermax
    % particle is on the right outside the buffer
    p_turn = k1*p0;
    %
elseif c_local > c_idealmax && c_local < c_buffermax && up_new*dc > 0
    % particle is between right edge and right buffer
    % wrong direction --> probability of turning is increased
    p_turn = k1*p0;
    %
elseif c_local < c_buffermin
    % particle is on the left outside the buffer
    p_turn = k2*p0;
    %
elseif c_local > c_buffermin && c_local < c_idealmin && up_new*dc < 0
    % particle is between left buffer and right edge
    % wrong direction --> probability of turning is increased
    p_turn = k2*p0;
    %
end

% applying the seed (different every time, but still reproducible)
rng('default');
rng(seed);
% computing the new direction, particle turns with a likelihood of p_turn
if rand() <= p_turn
    up_new = -up_new;
end

end