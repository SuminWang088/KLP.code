function [X, e, nEval] = blackboxseklpoints_greedy(nDim, fpdf, fmin, nIter,initial)
% [X, e, nEval] = blackboxseklpoints_greedy(nDim, fpdf, fmin, k, nIter) generates a point-
% set called adaptive Kullback-Leibler
% points (AKL points) using a one-point-at-a-time greedy for the black box function.
% 
%
% Input:
% nDim  - number of dimensions of the target density.
% fpdf  - handle to the target density function.
% fmin  - handle to a nDim-dimensional minimiser.
% nIter - length of the generated sequence of points.
% initial - the initial used to generate AKL points.
%
% Output:
% X     - nIter-by-nDim matrix of generated points.
% e     -the target function used to minimize at each iteration.
% nEval - number of density evaluations at each iteration.
%
%
    n0=size(initial,1);
    X = [initial;zeros((nIter-n0), nDim)];
    y = zeros(nIter, 1);
    e = zeros(nIter, 1);
    nEval = zeros(nIter, 1);

   
    
    for n = (n0+1):nIter
        f = @(XNew)efe(XNew, fpdf, X, y, n);
        [X(n, :), y(n), e(n), nEval(n)] = fmin(f, X(1:(n - 1), :));
        fprintf('n = %d\n', n);
    end
end

function [e, yNew] = efe(XNew, fpdf, X, y, n)
    [nNew, nDim] = size(XNew);
    S=X(1:(n - 1), :);
    Y=zeros(n-1,1);
    for i=1:(n-1)
        Y(i)=fpdf(S(i,:));
    end
    f0=@(XNew)dace(XNew,S,Y);
    x0=[0.2,0.2];
    lb = [0, 0];
ub = [1, 1];
    opt = optimset('tolfun', 1e-3, 'tolx', 1e-3, 'display', 'off');
    [~, f0min, ~, ~] = fminsearchcon( ...
            f0, x0, lb, ub, [], [], [], opt);
    yNew=dace(XNew,S,Y)-f0min;
    
    
    exist=X(1:(n - 1), :);
matrix=[XNew;exist];
d=pdist2(matrix,matrix);
[m,t]=size(d);
dnew=[];
for i=1:m
    dnew=[dnew;d(i,[1:i-1 i+1:t])];

end

de=pdist2(exist,exist);
[s,t]=size(de);
dd=[];
for i=1:s
    dd=[dd;de(i,[1:i-1 i+1:t])];

end
%%%%% The bandwidth h_n.
h_n=min(min(dd))/(n^(1/(nDim+4)));
dnew=exp(-dnew/h_n);

d1=sum(dnew(1,:));
he=log(d1);

ce=log(yNew);
c2=log(n*(h_n)^nDim);
kl=he-ce-c2;
e=kl;
    
end




     
