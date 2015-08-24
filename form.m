function mpfp = form(g,m,s,R,dists)

S = diag(s);

meq = zeros(length(m),1);
Seq = zeros(size(S));

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
C = Seq*R*Seq;
L = chol(C,'lower');
iL = inv(L);
Z = @(x) iL*(x-meq);
xz = @(z) L*z+meq;

dz = 1e-3;

zn = [0 0 0]';
crit = 1;
k=0;
while crit>1e-8
	zo = zn;
	X = xz(zo);
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

	grad = zeros(length(m),1);
	for i=1:length(m)
		zod = zo;
		zod(i) = zo(i)-dz;
		grad(i) = ( g(xz(zo))-g(xz(zod)) )/dz;
	end
	zn = grad * (grad'*zo-g(xz(zo)))/(grad'*grad);
	k=k+1;
	crit = sqrt((zn-zo)'*(zn-zo));
end

mpfp.k = k;
mpfp.x = xz(zn); % most probable failure point
mpfp.beta = sqrt(zn'*zn);
mpfp.Pf = normcdf(-mpfp.beta);

