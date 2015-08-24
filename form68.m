% exemplo 6.8 ang

g = @(x) x(1)*x(2)-x(3);                    % função de estado limite

M = [40 50 1000]';                          % vetor de médias
CV = [0.125 0.05 0.2]';                     % vetor de coeficientes de variação
S = M.*CV;                                  % vetor de desvios padrão
r = eye(3);                                 % matriz de correlação
dists = {'lognormal' 'lognormal' 'gumbel'}; % distribuições


mpfp = form(g,M,S,r,dists)                  % resultados
