% compute equivalent normal mean and deviation for
% lognormal distribution
function [mueq sigmaeq] = eqLN(x,mu,s)

zeta = sqrt(log(1+(s/mu)^2));
lambda = log(mu)-zeta^2/2;

sigmaeq = normpdf(norminv(logncdf(x,lambda,zeta)))/lognpdf(x,lambda,zeta);
mueq = x-sigmaeq*norminv(logncdf(x,lambda,zeta));
