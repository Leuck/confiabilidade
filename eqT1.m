% compute equivalent normal mean and deviation for
% type 1 extreme distribution
function [mueq sigmaeq] = eqT1(x,mu,s)

zeta = sqrt(log(1+(s/mu)^2));
lambda = log(mu)-zeta^2/2;

alpha = pi/(s*sqrt(6));
u = mu-0.5772/alpha;

t1cdf = @(x) exp(-exp(-alpha*(x-u)));
t1pdf = @(x) alpha*exp(-alpha*(x-u)-exp(-alpha*(x-u)));

sigmaeq = normpdf(norminv(t1cdf(x)))/t1pdf(x);
mueq = x-sigmaeq*norminv(t1cdf(x));
