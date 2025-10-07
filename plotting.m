function plotting(data,x0,channelLocations,name)

    steadyState = data(end,:);
    distances = channelLocations;

    p = polyfit(distances,steadyState,1);

    x1 = linspace(0,channelLocations(end));
    y1 = polyval(p,x1);
    T0 = polyval(p,x0);

    figure()
    scatter(distances,steadyState);
    hold on;
    plot(x1,y1);
    xlabel('Location (m)')
    ylabel('Temperature (deg C)')
    legend('Raw Data','Best Fit')
    title(name)

    fprintf("For the %s case, T0 = %f \n",name,T0);
    fprintf("For the %s case, H_exp = %f \n",name,p(1));


end