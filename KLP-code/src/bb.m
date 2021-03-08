function f=bb(x)
% The two-dimensional black box function described by Fasshauer (2007).
  f=3*exp(-(9*x(1)-2)^2/4-(9*x(2)-2)^2/4)/4+3*exp(-(9*x(1)+1)^2/49-(9*x(2)+1)^2/10)/4+exp(-(9*x(1)-7)^2/4-(9*x(2)-3)^2/4)/2-exp(-(9*x(1)-4)^2-(9*x(2)-7)^2)/5;
end

