function [fitresult, gof] = LinearFit(xData, yData)

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData(xData, yData);

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );





