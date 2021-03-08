function [X, e, nEval] = aklpoints_greedy(nDim, fpdf, fmin, nIter,initial)
% [X, e, nEval] = AKL_greedy(nDim, fpdf, fmin,  nIter,initial) generates a point-
% set called adaptive Kullback-Leibler points (AKL) using a one-point-at-a-time greedy
% algorithm.
%
% Input:
% nDim  - number of dimensions of the target density.
% fpdf  - handle to the target density function.
% fmin  - handle to a nDim-dimensional minimiser.
% nIter - the number of AKL points.
% initial- the initial set used to generate AKL points.
% Output:
% X     - nIter-by-nDim matrix of generated AKL points.
% e     - the Kullback-Leibler divergence brtween the kernel density of AKL
%%%% points and the Gaussian process model at each iteration.
% nEval - number of density evaluations at each iteration.
%
%
    n0=size(initial,1);
    X = [initial;zeros((nIter-n0), nDim)];
    y = zeros(nIter, 1);
    e = zeros(nIter, 1);
    nEval = zeros(nIter, 1);

   
    
    for n = (n0+1):nIter
        f = @(XNew)efe(XNew, fpdf, X,  n);
        [X(n, :), y(n), e(n), nEval(n)] = fmin(f, X(1:(n - 1), :));
        fprintf('n = %d\n', n);
    end
end

function [e, yNew] = efe(XNew, fpdf, X, n)
    [~, nDim] = size(XNew);
    %%Constructing Gaussian process model based on the first n-1 points;
    S=X(1:(n - 1), :);
    Y=zeros(n-1,1);
    for i=1:(n-1)
        Y(i)=log(fpdf(S(i,:)));
    end
    yNew=dace(XNew,S,Y);%%%Prediction function;
    %%%%%%%%%%%%%%%%%%
    
    %%%%%% the matrix of current point;
    exist=X(1:(n - 1), :);
    %%%%%%%
matrix=[XNew;exist];
d=pdist2(matrix,matrix);
[m,t]=size(d);
dnew1=[];
for i=1:m
    dnew1=[dnew1;d(i,[1:i-1 i+1:t])];

end
de=pdist2(exist,exist);
[s,t]=size(de);
dd=[];
for i=1:s
    dd=[dd;de(i,[1:i-1 i+1:t])];

end
%%%%% The bandwidth h_n
h_n=min(min(dd))./(n^(1/(nDim+4)));
%%%%%%%%

dnew=exp(-dnew1(1,:)/h_n);

    d1=sum(dnew(1,:)) ;

he=sum(log(d1));


pre=yNew;%%%% Prediction of the log-density.

kl=he-pre-log(n*(h_n)^nDim);
%%%% e is the Kullback-Leibler divergence brtween the kernel density of AKL
%%%% points and the Gaussian process model;

e=kl;
    
    
end



     
