function lines(c0,x0)

% this function only serves the purpose to visualize the ideal
% concentration and the buffer region for aerokinesis. It is only called
% inside "harryplotter.m" and has no further relevance...

% parameters for the lines. Need to be changed according to the ones in
% Aerotaxis to prevent visual misinterpretations!!
c_idealmin = 0.4*c0;
c_idealmax = 0.6*c0;
c_buffermin = 0.3*c0;
c_buffermax = 0.8*c0;

% here, we create the four vertical lines. That´s it.
hold on
xline(x0*log(c0/c_idealmax),'r')
xline(x0*log(c0/c_idealmin),'r')
xline(x0*log(c0/c_buffermax))
xline(x0*log(c0/c_buffermin))
hold off

end