% montecarlo(g,m,s,R,n,dists)
%
% Entradas:
%   g: handle para a função de estado limite
%   m: vetor de médias
%   s: vetor de desvios padrão
%   R: matriz de correlação
%   n: número de pontos aleatórios para avaliar a f.e.l.
%   dists: nome das distribuições
% Saídas:
%   O: estrutura contendo a probabilidade de falha (Pf), o coeficiente de
%      variação (CV) e o beta (beta).
function O = montecarlo(g,m,s,R,n,dists)

S = diag(s);         % matriz S
C = S*R*S;           % matriz de covariância
L = chol(R,'lower'); % fatoração da matriz de correlação

% n conjuntos de variáveis normais não correlacionadas
Z = normrnd(0,1,[length(m) n]);
Zc = L*Z;           % imposição de correlação
Uc = normcdf(Zc);   % variáveis correlacionadas distribuidas uniformemente

% variáveis correlacionadas nas distribuições desejadas
X = zeros(size(Z));
for i=1:length(m)
	if strcmp(dists{i},'normal')
		X(i,:)= norminv(Uc(i,:),m(i),s(i));
	elseif strcmp(dists{i},'lognormal')
		zeta = sqrt(log(1+(s(i)/m(i))^2));
		lambda = log(m(i))-zeta^2/2;
		X(i,:)= logninv(Uc(i,:),lambda,zeta);
	elseif strcmp(dists{i},'gumbel')
		alpha = pi/(s(i)*sqrt(6));
		u = m(i)-0.5772/alpha;
		X(i,:)= (-log(-log(Uc(i,:)))+alpha*u)/alpha;
	end
end

Pf = sum(g(X)<0)/n;         % probabilidade de falha
O.Pf = Pf;
O.CV = sqrt((1-Pf)/(n*Pf)); % coeficiente de variação
O.beta = norminv(1-Pf);     % beta

