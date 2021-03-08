function d = lhd(n, p)
d=zeros(n,p);
for i=1:p
    d(:,i)=randperm(n)';
end
    d=(d-rand(n,p))/n;
end
