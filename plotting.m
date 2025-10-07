function plotting(data,x0,channelLocations)

    steadyState = data(end,:);
    distances = [x0,channelLocations];

    figure()
    scatter(distances,steadyState);


end