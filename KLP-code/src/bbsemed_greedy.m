function [X, e, nEval] = bbsemed_greedy(nDim, fpdf, fmin, k, nIter,initial)
% [X, e, nEval] = bbsemed_greedy(nDim, fpdf, fmin, k, nIter,initial) generates a point-
% set called sequential minimum energy design (SMED) using a
% one-point-at-a-time greedy for the black box function
% algorithm described by Joseph et al (2015).
%
% Input:
% nDim  - number of dimensions of the target density.
% fpdf  - handle to the target density function.
% fmin  - handle to a nDim-dimensional minimiser.
% k     - power parameter of the generalised energy criterion.
% nIter - length of the generated sequence of points.
% initial - the initial used to generate SMED points.
%
% Output:
% X     - nIter-by-nDim matrix of generated points.
% e     - minimised energy at each iteration.
% nEval - number of density evaluations at each iteration.
%
%

   n0=size(initial,1);
    X = [initial;zeros((nIter-n0), nDim)];
    y = zeros(nIter, 1);
     for i=1:n0
        y(i)=fpdf(X(i,:)) .^ (k ./ (2 .* nDim));
    end
    e = zeros(nIter, 1);
    nEval = zeros(nIter, 1);

  
    
    for n = (n0+1):nIter
        f = @(XNew)efe(XNew, fpdf, k, X, y, n);
        [X(n, :), y(n), e(n), nEval(n)] = fmin(f, X(1:(n - 1), :));
        fprintf('n = %d\n', n);
    end
end



function [e, yNew] = efe(XNew, fpdf, k, X, y, n)
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
        
    fpdfy=dace(XNew,S,Y)-f0min;
    yNew = fpdfy .^ (k ./ (2 .* nDim));
    A = repmat(XNew, n - 1, 1);
    B = repelem(X(1:(n - 1), :), nNew, 1);
    yb = repelem(y(1:(n - 1)), nNew, 1);
    t = 1 ./ (yb .* fd(A, B) .^ k);
    e = 1 ./ yNew .* sum(reshape(t, nNew, []), 2);
end

function [q, yNew] = fq(XNew, fpdf)
    nDim = size(XNew, 2);
    yNew = fpdf(XNew) .^ (1 ./ (2 .* nDim));
    q = 1 ./ yNew;
end

function d = fd(A, B)
    d = sqrt(sum((A - B) .^ 2, 2));
end

