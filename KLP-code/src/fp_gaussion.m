function fx = fp_gaussion(p,x)
u=zeros(1,p);
sigma=zeros(p,p);
for i=1:p
    sigma(i,i)=1;
end
fx=exp(-(x-u)*inv(sigma)*(x-u)')/((2*pi)^(p/2)*sqrt(det(sigma)));
end

