function [X, e, nEval] = klpointsnr05_greedy(nDim, fpdf, fmin, nIter)
% [X, e, nEval] = klpointsnr05_greedy(nDim, fpdf, fmin, nIter) generates a point-
% set called Kullback-Leibler points (KL) using a one-point-at-a-time greedy
% for the  bivariate normal distribution N(0,R^2), where R is a 2*2
% matrix with the R(i,j)=0.5^{|i-j|};
%
% Input:
% nDim  - number of dimensions of the target density.
% fpdf  - handle to the target density function.
% fmin  - handle to a nDim-dimensional minimiser.

%  nIter -the number of KL points.
%
% Output:
% X     - nIter-by-nDim matrix of generated points.
% e     - the target function used to minimize at each iteration.
% nEval - number of density evaluations at each iteration.
%
%

    X = zeros(nIter, nDim);
    y = zeros(nIter, 1);
    e = zeros(nIter, 1);
    nEval = zeros(nIter, 1);

    % Generate x_1
    f = @(XNew)fq(XNew, fpdf);
    [X(1, :), y(1), e(1), nEval(1)] = fmin(f, double.empty(0, nDim));
    fprintf('n = 1\n');

    % Generate the rest
    for n = 2:nIter
        f = @(XNew)fe(XNew, fpdf, X, y, n);
        [X(n, :), y(n), e(n), nEval(n)] = fmin(f, X(1:(n - 1), :));
        fprintf('n = %d\n', n);
    end
end

function [e, yNew] = fe(XNew, fpdf, X, y, n)
    [nNew, nDim] = size(XNew);
    yNew = fpdf(XNew);
    exist=X(1:(n - 1), :);
matrix=[XNew;exist];
d=pdist2(matrix,matrix);
[m,t]=size(d);
dnew=[];
for i=1:m
    dnew=[dnew;d(i,[1:i-1 i+1:t])];

end
c1=n^(-1/(nDim+4));
%%%%% The bandwidth h_n.
h_n=0.8*c1;
dnew=exp(-dnew/h_n);

d1=sum(dnew(1,:));
he=log(d1);

ce=log(yNew);
c2=log(n*(h_n)^nDim);
kl=he-ce-log(c2);
e=kl;
    
    
end


function [q, yNew] = fq(XNew, fpdf)
    yNew = fpdf(XNew);
    q = 1 ./ yNew;
end
