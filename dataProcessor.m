function cleanData = dataProcessor(data)

    cleanData = data(~any(isnan(data), 2), :);
    cleanData(:,1) = [];


end