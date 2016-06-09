function [y_m, R, P, b, y_int] = LinearRegression(x, y, ifPlot)
% fits x and y to a linear model, and returns the R and P values, along
% with y_m, model y values
%
%  [y_m, R, P, b, y_int] = LinearRegression(x, y)

degree=2; % linear regression

if length(x)~=length(y)
    error('LinearRegression needs x and y with same lengths');
end

if length(x)==1
    disp('Only one data point to LinearRegression, so 1s returned');
    y_m=1;
    R=0;
    P=0;
    return;
end

tf = isnan(x) | isnan(y); % remove NaNs
x(tf) = [];
y(tf) = [];

lincoeff = polyfit(x, y, 1);

y_m = polyval(lincoeff, x, degree);

model_v_observed(:,1) = y;
model_v_observed(:,2) = y_m;

[r,p]=corrcoef(model_v_observed); %to determine p values

R = r(1,2);
P = p(1,2);

b = lincoeff(1);

y_int = lincoeff(2);

if ~exist('ifPlot','var')
    ifPlot = 0;
end

if ifPlot==1
    
    plot(x,y,'.')
    hold on
    xv = [min(x) max(x)];
    plot(xv,xv*b+y_int);
    title(['R^2=' num2str(R), ' p=' num2str(P)])
end