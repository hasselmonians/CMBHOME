function [y_model,fit_params,rsquare,goodness] = SaturatingRegression(x, y)
% Does a natural growth fit and returns the model expected values,
% parameters of best fit, rsquare value, and other goodness object.
%
% [fit,params,rsquare,goodness] = CMBHOME.Utils.LogRegression(x,y)

if length(x)~=length(y)
    error('LogRegression needs x and y with same lengths');
end

if length(x)==1
    disp('Only one data point to LogRegression, so 1s returned');
    y_m=1;
    R=0;
    P=0;
    return;
end

x = x(:);1
y = y(:);

tf = isnan(x) | isnan(y); % remove NaNs
x(tf) = [];
y(tf) = [];

f = fittype('d - a*exp(-b*x)','independent','x','coefficients',{'d' 'a' 'b'});
options = fitoptions(f);
options.Lower = [0 0 0];
options.Upper = [Inf,Inf,Inf];

[fit_params, goodness] = fit(x,y,f,options);
rsquare = goodness.rsquare;

y_model = fit_params.d - fit_params.a*exp(-fit_params.b*x);