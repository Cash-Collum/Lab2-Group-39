function plotting(data,x0,channelLocations,name)

    steadyState = data(end,:);
    distances = [x0,channelLocations];

    p = polyfit(distances,steadyState,1);

    x1 = linspace(0,channelLocations(end));
    y1 = polyval(p,x1);

    figure()
    scatter(distances,steadyState);
    hold on;
    plot(x1,y1);
    xlabel('Location (m)')
    ylabel('Temperature (deg C)')
    legend('Raw Data','Best Fit')
    title(name)

    fprintf("For the %s case, T0 = %f \n",name,steadyState(1));


end