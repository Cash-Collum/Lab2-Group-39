function plotting(steadyState,initState,distances,x1,y1,x2,y2,analytical,linear,raw,name)

    figure()
    hold on;
    scatter(distances,steadyState);
    scatter(analytical.xTherm,raw);
    plot(x1,y1);
    plot(analytical.x,linear);
    xlabel('Location (m)')
    ylabel('Temperature (deg C)')
    legend('EXP Raw Data','AN Raw Data','EXP Linear','AN Linear')
    title("Steady State: " + name)
    hold off;

    saveas(gcf, "Steady_State " +name + ".png");


    figure()
    scatter(distances,initState);
    hold on;
    plot(x2,y2);
    xlabel('Location (m)')
    ylabel('Temperature (deg C)')
    legend('Raw Data','Best Fit')
    title("Initial State: " + name)

    saveas(gcf, "Initial_State " +name + ".png");


end