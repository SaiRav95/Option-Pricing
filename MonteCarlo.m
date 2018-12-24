%% Input parameters
S = 100;    %stock price at time 0
K = 90;     %strike price
T = 1;      %maturity in years
s = 0.4;    %annual volatility
r = 0.06;   %continuously compounded annual interest rate
n = 365;     %number of periods
k = 10000; %number of simulation paths
dt = T/n;

%generating sample paths
nu = (r - 0.5*s^2)*dt;
si = s*sqrt(dt);
randomMatrix = randn(k, n);
increments = nu+si*randomMatrix;

logpaths = cumsum([log(S)*ones(k, 1), increments], 2);
pathMatrix = exp(logpaths);

for line = 1:100
    plot(pathMatrix(line,:))
    hold on
end

CashFlow = exp(-r*T) * max(0, pathMatrix(:,end)-K);
price = mean(CashFlow);
se = std(CashFlow)/sqrt(k);
blprice = blsprice(S, K, r, T, s);

fprintf('Option Price: %.3f\n', price);
fprintf('Standard Error: %.3f\n', se);
fprintf('Black-Scholes Price: %.3f\n', blprice);