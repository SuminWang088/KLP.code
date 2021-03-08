function y=dace(x,S,Y)
[nrow,ndim]=size(x);
  theta =10*ones(1,ndim); lob = 1e-1*ones(1,ndim); upb = 20*ones(1,ndim);
  [dmodel, perf] = dacefit(S, Y, @regpoly0, @corrgauss, theta, lob, upb);
  [YX,MSE] = predictor(x, dmodel);
  y=YX;     
end

