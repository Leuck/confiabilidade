% exemplo 6.9 ang
% metodo form

g = @(x) x(1)*x(2)-x(3);

M = [40 50 1000]';
CV = [0.125 0.05 0.2]';
S = M.*CV;
r = [1 .4 0; .4 1 0; 0 0 1];
dists = {'normal' 'normal' 'normal'};

mpfp = form(g,M,S,r,dists)
