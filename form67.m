% exemplo 6.7 ang
% método FORM

g = @(x) x(1)*x(2)-x(3);              % função de estado limite

M = [40 50 1000]';                    % vetor de médias
CV = [0.125 0.05 0.2]';               % vetor de coeficientes de variação
S = M.*CV;                            % vetor de desvios padrão
r = eye(3);                           % matriz de correlação
dists = {'normal' 'normal' 'normal'}; % distribuições

mpfp = form(g,M,S,r,dists)            % resultados
