function [mdl] = HyperbolicFit(xData, yData, startPoints)

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData(xData, yData);

modelfun = @(b,x) (b(1)*x(:, 1))./(x(:, 1) + b(2));
tbl = table(xData,yData);
mdl = fitnlm(tbl,modelfun,startPoints);

% Set up fittype and options.
% ft = fittype( '(a*x)/(x+b)', 'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Algorithm = 'Levenberg-Marquardt';
% opts.Display = 'Off';
% opts.StartPoint = startPoints;
% 
% % Fit model to data.
% [fitresult, gof, output] = fit( xData, yData, ft, opts );




