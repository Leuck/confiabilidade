function Pf = montecarlo(g,m,s,R,n,dists)
%n=1e6;

S = diag(s);
C = S*R*S;
L = chol(R,'lower');

Z = normrnd(0,1,[length(m) n]);
Zc = L*Z;
Uc = normcdf(Zc);
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

Pf = sum(g(X)<0)/n
CVPf = sqrt((1-Pf)/(n*Pf))
beta = norminv(1-Pf)
