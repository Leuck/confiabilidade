% exemplo 6.9 ang

% função de estado limite
g = @(x) x(1,:).*x(2,:)-x(3,:);

M = [40 50 1000]';
CV = [0.125 0.05 0.2]';
S = M.*CV;
r = [1 .4 0; .4 1 0; 0 0 1];
dists = {'normal' 'normal' 'normal'};

n = 4e6;

O = montecarlo(g,M,S,r,n,dists)
