
e=zeros(1000,1);
for j=1:1000
    nPart=8000;
X = gaussmix_rnd(Mu, C, w, nPart);
n=nPart;
mcg=zeros(n,1);
for i=1:n
    mcg(i)=gs(X_kl(i,:));
end
mean=abs(sum(mcg)/n-0.1148);
e(j)=mean;
end
log(sum(e)/1000)



%%%%%计算kl点的积分误差%%%%%
e=zeros(50,1);
for j=1:50
nPart=100;
[X_kl, ~, nEval_11] = klpointsk_greedy(d, fp, fmin_1, nPart);
n=nPart;
mcg=zeros(n,1);
for i=1:n
    mcg(i)=gosc(X_kl(i,:));
end
mean=abs(sum(mcg)/n);
e(j)=mean;
end
log(sum(e)/50)


%%%%%计算med点的积分误差%%%%%
e=zeros(50,1);
for j=1:50
nPart=50;
[X_kl, ~, nEval_11] = klpointsk_greedy(d, fp, fmin_1, nPart);
n=nPart;
mcg=zeros(n,1);
for i=1:n
    mcg(i)=g(X_kl(i,:));
end
mean=sum(mcg)/n;
e(j)=log(mean);
end
sum(e)/50