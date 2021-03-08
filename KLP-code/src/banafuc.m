function f=banafuc(x)
% The banana-shaped density function discribed by Haario et al. (1999).
  l1 = -20;
  u1=20;
  l2 = -10;
  u2 = 5;
  x1 = l1 + (u1 - l1)*x(1);
  x2 =l2 + (u2 - l2)*x(2);
  f=exp(-0.5*(x1 ^2 / 100 + (x2+ 0.03*x1^2 -3)^2));
end


