function harryplotter(type,wall,pd_max,x,v,vmax,vrange,u,sigma_pot,npart,c0,x0,v0,c_plat)

% this is the plotting function, which serves for the Programs General (2),
% Chemokinesis (3) and Aerotaxis (3). The initial plot (1) is always the
% same. The integer type is used to switch between those cases.

% in case there is no concentration gradient, the respective parameters are
% all set to zero.
if nargin < 11
    c0 = 0;
    x0 = 0;
end

% in case there is no velocity gradient, the respective parameters are
% all set to zero.
if nargin < 13
    v0 = 0;
    c_plat = 0;
end

% setting the bin width for the histogram
BinW=1/2*wall/floor(sqrt(npart));

% checking whether to make initial or final plot
switch type
    case 1
        figure(1)
    case 2
        figure(2)
    case 3
        figure(2)
    case 4
        figure(2)
end

% 
[nhist, edges] = histcounts(x,'BinWidth',BinW);
bins = zeros(1,length(edges)-1);
for ib = 1:(length(edges)-1)
    bins(ib) = edges(ib) + (edges(ib+1) - edges(ib))/2;
end
yhist = nhist/(BinW*npart);

% 
s = spline(bins,yhist);

% make the first subplot
subplot(3,1,1);

%% artefacts
%histogram(x,'BinWidth',1/2*BinW,'Normalization','probability')
%fit_pd = fitdist(x','normal');
%pd = normpdf(x,fit_pd.mu,fit_pd.sigma);
%plot(x,pd,'LineWidth',3)
%%

% plot the histogram of the particle distribution into the upper subplot, 
% and the spline on top of it
histogram(x,'BinWidth',BinW,'Normalization','pdf','FaceColor','b')
hold on
plot(bins,yhist,bins,ppval(s,bins),'color','b','linestyle','-','LineWidth',2)
hold off

% assign the right title (initial or real time)
switch type
    case 1 % initial
        title("initial particle distribution")
    case 2 % general
        title("real time particle distribution")
    case 3 % kinesis
        title("real time particle distribution")
    case 4 % taxis
        title("real time particle distribution")
        % plot vertical lines for visualizing the ideal concentration
        % (red) and the buffer region (black)
        lines(c0,x0)
end

% adjusting labels and x/y-ranges of the plot
xlabel("particle positions")
ylabel("normalized count")
xlim([0 wall])
ylim([0 pd_max])
axis('normal')

subplot(3,1,2)

% 
interval = 0:wall;
switch type
    case 1 % initial
        plot(x,v)
        title("initial velocity distribution")
        ylim([0 vrange])

    case 2 % general
        % plot the velocities of all the particles and the mean velocity.
        % velocities should be distributed around it inside the range.
        plot(x,v)
        hold on
        plot(interval,0.5*vrange,'r')
        hold off
        title("real time velocity distribution")
        ylim([0 vmax])

    case 3 % kinesis
        % plot the velocities of all the particles
        plot(x,v)

        % calculating all the local concentrations and putting them in one
        % array c_global. This is necessary to calculate the velocity
        % gradient in the next step.
        c_global = local_c(wall-interval,c0,x0);

        % calculating the local velocities for each local concentration
        % (from c_global). Then assigning it to the global velocity array 
        % v_global for plotting.
        %%This is indeed quite complicated. There might be an easier way...
        v_global = zeros(1,length(c_global));
        ii = 1;
        for c=c_global
            v_global(ii) = local_v(0,c,c_plat,v0);
            ii=ii+1;
        end
        
        % plotting the velocity distribution v_global on top of the actual
        % velocities.
        hold on
        plot(interval,v_global,'r')
        legend("","velocity gradient")   
        hold off
        % setting the right title and y-range
        title("real time velocity distribution")
        ylim([0 1.7*vmax])

    case 4 % taxis
        % separating between left and right direction to assign the right
        % marker (< or >, blue or red) in the plots
        nright = nnz(u+1);
        nleft = nnz(u-1);
        right = zeros(2,nright);
        left = zeros(2,nleft);
        index = 1:npart;
        for pp = index
            if u(pp) == 1
                right(:,pp) = [x(pp);v(pp)];
            else
                left(:,pp) = [x(pp);v(pp)];
            end
        end
        % plot all particles with left (blue) and right (red) direction,
        % adjust the x/y-range
        plot(nonzeros(right(1,:)),nonzeros(right(2,:)),'r>')
        hold on
        plot(nonzeros(left(1,:)),nonzeros(left(2,:)),'b<')
        hold off
        title("real time velocity distribution")
        xlim([0 wall])
        ylim([0 vmax])
        % plot vertical lines for visualizing the ideal concentration
        % (red) and the buffer region (black)
        lines(c0,x0)
end

% plot either particle positions (1: initial or 2: general)
if type <= 2
    xlabel("particle positions")
    ylabel("velocity of particles")
    xlim([0 wall])
    % ylim([0 vmax])
    axis('normal')

    % create the bottom subplot (1D/axes plot) for the particle positions
    s = subplot(3,1,3);
    set(s,'visible','off'); %make the additional plot behind invisible
    hAxes = axes('NextPlot','add',... % Add subsequent plots to the axes,
        'DataAspectRatio',[1 1 1],... % match the scaling of each axis,
        'XLim',[0 wall],... % set the x axis limit,
        'YLim',[-eps eps]); % set the y axis limit (tiny!)
    set(hAxes,'position',[0.13,0.15,0.77,0.15]);
    plot(x,0,'.','MarkerSize',25*sigma_pot); %scale it properly!
    %axis('equal')

% or concentration profile (3: chemokinesis or 4: aerotaxis)
else
    subplot(3,1,3)
    xlabel("particle positions")
    ylabel("concentration of oxygen")
    xlim([0 wall])
    variable = 0:wall;
    plot(variable,c0*exp(-variable/x0))
end

% ... and finally, set the right title. Now, all the subplots are ready!
switch type
    case 1 % initial
        titl = title("initial particle positions");
        set(titl,'position',[wall/2, 1000*eps]);
    case 2 % general
        titl = title("real time particle positions");
        set(titl,'position',[wall/2, 1000*eps]);
    case 3 % kinesis
        title("concentration profile");
    case 4 % taxis
        title("concentration profile");
        % plot vertical lines for visualizing the ideal concentration
        % (red) and the buffer region (black)
        lines(c0,x0)
end

