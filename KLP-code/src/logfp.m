function f=fp(x,d)
sigma=eye(d);
mu=zeros(1,d);
f =(2*pi)^(-1/(2*d))*det(sigma)^(-1/2)*exp(-(x-mu)*inv(sigma)*(x-mu)'/2);
  
end




