function output = getSaturation(x0,maxCurveValue,logisticGrowth,x)

% https://en.wikipedia.org/wiki/Logistic_function

%INPUT----
%- x0 (double, [-])                       the x value of the sigmoid's midpoint;
%- maxCurveValue (double, [-])            the curve's maximum value;
%- logisticGrowth (double, [-])           the logistic growth rate or steepness of the curve.
%- x (double, [-])                        independent variable

%OUTPUT----
%- output (double,[-])                    Returns function value

output = maxCurveValue ./ ( 1 + exp(-logisticGrowth .* ( x - x0 )) );

end
