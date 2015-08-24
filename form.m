% form(g,m,s,R,dists)
%
% Entradas:
%   g: handle para a função de estado limite
%   m: vetor de médias
%   s: vetor de desvios padrão
%   R: matriz de correlação
%   dists: nome das distribuições
% Saídas:
%   mpfp: estrutura contendo o numero de iterações (k), o ponto de falha mais
%         provável (x), o beta (beta) e a probabilidade de falha (Pf).
function mpfp = form(g,m,s,R,dists)

S = diag(s); % matriz S

meq = zeros(length(m),1);
Seq = zeros(size(S));

% para cada variavel calcular media e desvio padrão normal equivalente
for i=1:length(m)
	if strcmp(dists{i},'normal')
		Seq(i,i)=S(i,i);
		meq(i)=m(i);
	elseif strcmp(dists{i},'lognormal')
		[meq(i) Seq(i,i)] = eqLN(m(i),m(i),S(i,i));
	elseif strcmp(dists{i},'gumbel')
		[meq(i) Seq(i,i)] = eqT1(m(i),m(i),S(i,i));
	end
end
C = Seq*R*Seq;        % Matriz de covariancia
L = chol(C,'lower');  % Matriz L
iL = inv(L);          % inversa de L
Z = @(x) iL*(x-meq);  % função Z(X)
xz = @(z) L*z+meq;    % função X(Z)

dz = 1e-3;            % dz para calculo dos gradientes (diferenças finitas)

zn = [0 0 0]';        % estimativa inicial
crit = 1;
k=0;

%%%%%%%%%%%%%%%%%%%%
%% LOOP PRINCIPAL %%
%%%%%%%%%%%%%%%%%%%%

while crit>1e-8 % critério de parada
	zo = zn;
	X = xz(zo);
	% para cada variavel calcular media e desvio padrão normal equivalente
	for i=1:length(m)
		if strcmp(dists{i},'normal')
			Seq(i,i)=S(i,i);
			meq(i)=m(i);
		elseif strcmp(dists{i},'lognormal')
			[meq(i) Seq(i,i)] = eqLN(X(i),m(i),S(i,i));
		elseif strcmp(dists{i},'gumbel')
			[meq(i) Seq(i,i)] = eqT1(X(i),m(i),S(i,i));
		end
	end
	C = Seq*R*Seq;
	L = chol(C,'lower');
	iL = inv(L);
	Z = @(x) iL*(x-meq);
	xz = @(z) L*z+meq;

	% calculo do gradiente
	grad = zeros(length(m),1);
	for i=1:length(m)
		zod = zo;
		zod(i) = zo(i)-dz;
		grad(i) = ( g(xz(zo))-g(xz(zod)) )/dz;
	end
	% novo vetor de variaveis normalizadas
	zn = grad * (grad'*zo-g(xz(zo)))/(grad'*grad);
	k=k+1;
	% mudança no vetor de variaveis normalizadas
	crit = sqrt((zn-zo)'*(zn-zo));
end

mpfp.k = k;                    % número de iterações
mpfp.x = xz(zn);               % most probable failure point
mpfp.beta = sqrt(zn'*zn);      % beta
mpfp.Pf = normcdf(-mpfp.beta); % probabilidade de falha

