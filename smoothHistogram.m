function smoothHistogram(data, numBins,color)
    % Create a histogram
    [counts, edges] = histcounts(data, numBins, 'Normalization', 'pdf');

    % Calculate bin centers
    binCenters = (edges(1:end-1) + edges(2:end)) / 2;
    
    
    % Use spline to create a smooth curve
    pp = spline(binCenters, counts);

    % Generate a finer set of x-values
    xFine = linspace(min(binCenters), max(binCenters), 100);

    % Evaluate the spline at the finer x-values
    ySmooth = ppval(pp, xFine);

    delta = 0.01;
    xFine = [min(xFine)-delta xFine max(xFine)+delta];
    ySmooth = [0 ySmooth 0];
    ySmooth(ySmooth < 0) = 0;
    
    % Plot the smooth histogram
    plot(xFine, ySmooth, 'LineWidth', 2,Color= color);
    %plot(binExtended,countsExtended,'LineWidth', 2);
end


% function smoothHistogram(data, numBins)
%     % Create a histogram
%     [counts, edges] = histcounts(data, numBins, 'Normalization', 'pdf');
% 
%     % Calculate bin centers
%     binCenters = (edges(1:end-1) + edges(2:end)) / 2;
%     
%     delta = 0.01;
%     % Extend the data and counts with control points to force the spline to zero
%     binExtended = [min(binCenters)-delta binCenters max(binCenters)+delta];
%     countsExtended = [0 counts 0];
% 
%     % Use spline to create a smooth curve
%     pp = spline(binExtended, countsExtended);
% 
%     % Generate a finer set of x-values
%     xFine = linspace(min(binCenters), max(binCenters), 100);
% 
%     % Evaluate the spline at the finer x-values
%     ySmooth = ppval(pp, xFine);
% 
% 
% 
%     % Plot the smooth histogram
%     %plot(xFine, ySmooth, 'LineWidth', 2);
%     plot(binExtended,countsExtended,'LineWidth', 2);
% end
% 


