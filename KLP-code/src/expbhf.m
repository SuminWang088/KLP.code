function f=expbhf(x)
  l1 = -5;
  u1=10;
  l2 = 0;
  u2 = 15;
  x1 = l1 + (u1 - l1)*x(1);
  x2 =l2 + (u2 - l2)*x(2);
  f1=(x2-5.1*x1^2/(4*pi^2)+5*x1/pi-6)^2+10*((1-1/(8*pi))*cos(x1)+1)+5*x1;
  f=exp(f1);
end
