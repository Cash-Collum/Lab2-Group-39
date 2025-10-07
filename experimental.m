function [steadyState,initState,x1,y1,x2,y2] = experimental(data,x0,channelLocations,name)

    steadyState = data(end,:);
    initState = data(1,:);

    p = polyfit(channelLocations,steadyState,1);
    p2 = polyfit(channelLocations,initState,1);

    x1 = linspace(0,channelLocations(end));
    y1 = polyval(p,x1);
    T0 = polyval(p,x0);

    x2 = linspace(0,channelLocations(end));
    y2 = polyval(p2,x2);

    fprintf("For the %s case, T0 = %f \n",name,T0);
    fprintf("For the %s case, H_exp = %f \n",name,p(1));

    fprintf("For the %s case, M_exp = %f \n",name,p2(1));


end