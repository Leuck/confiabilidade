% exemplo 6.10 ang

g = @(x) x(1,:).*x(2,:)-x(3,:);

M = [40 50 1000]';
CV = [0.125 0.05 0.2]';
S = M.*CV;
%r = eye(3);
r = [1 .4 0; .4 1 0; 0 0 1];

%dists = {'normal' 'normal' 'normal'};
dists = {'lognormal' 'lognormal' 'gumbel'};

n = 2e6;

X = montecarlo(g,M,S,r,n,dists);
