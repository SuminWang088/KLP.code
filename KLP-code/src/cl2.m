function cl=cl2(x)
[n,p]=size(x);
pa1=0;
for k=1:n
    for i=1:p
        p1=1+0.5*abs(x(k,i)-0.5)-0.5*abs(x(k,i)-0.5)^2;
        p1=p1*p1;
    end
     pa1=pa1+p1;
end
part1=2*pa1/n;

pa2=0;
for k=1:n
    for j=1:n
        for i=1:p
            p2=1+0.5*abs(x(k,i)-0.5)+0.5*abs(x(j,i)-0.5)-0.5*abs(x(k,i)-x(j,i));
           p2=p2*p2;
        end
        pa2=pa2+p2;
    end
end
part2=pa2/n^2;
cl=(13/12)^p-part1+part2;
end


