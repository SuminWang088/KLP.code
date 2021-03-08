
fp = @(X)bhf(X);
fpp=@(X)expbhf(X);
fu = @(X)fscr_gaussmix(X, Mu, C, w);

% Symbolic variables
a = sym('a', [1, 2], 'real');
b = sym('b', [1, 2], 'real');

% Optimisers
lb = [0, 0];
ub = [1, 1];
mu0 = [0, 0];
Sigma0 = diag([25, 25]);
nStep = 100;
fmin_1 = @(f, X)fmin_adares(f, X, 2, mu0, Sigma0, 1, 30, lb, ub, [], []);
fmin_2 = @(f, X)fmin_adamc(f, X, 20, mu0, Sigma0, 1, 20, lb, ub, [], []);
fmin_3 = @(f, X)fmin_gsrtn(f, X, lb, ub, [nStep, nStep]);
fmin_4 = @(f, X)fmin_gs2d(f, lb, ub, [nStep, nStep]);
fmin_5 = @(f, X)fmin_ps(f, 2, 10, lb, ub);
fmin = fmin_2;

% Extensible point sets
nPart_1 = 8000;
nPart_2 = 500;
nPart_3 = 500;
nPart_4 = 500;

% (1) Monte Carlo
X = gaussmix_rnd(Mu, C, w, nPart_1);
nEval = ones(nPart_1, 1);
X_1 = cumset(X);
nEval_1 = log(cumsum(nEval));
nPart=109;
[X_med, ~, nEval_1] = med_greedy(2, fp, fmin_1, k_1, nPart);
[X_med, ~, nEval_1] = smed_greedy(2, fp, fmin_1, k_1, nPart);
[X_med, ~, nEval_1] = semed_greedy(2, fp, fmin_1, k_1, nPart,initial1);
[X_med, ~, nEval_1] = bbsemed_greedy(2, fp, fmin_1, k_1, nPart,initial);
[X_med, ~, nEval_1] = semed_greedyjh(2, fp, fmin_1, k_1, nPart,initial);

[X_kl, ~, nEval_10] = klpoints_greedy(2, fp, fmin_1, nPart);
[X_kl, ~, nEval_10] = klpointsk_greedy(2, fp, fmin_1, nPart);
[X_klnew, ~, nEval_11] = klpointsnew_greedy(2, fp, fmin_1, nPart);
[X_skl, ~, nEval_12] = sklpoints_greedy(2, fp, fmin_1, nPart);
[X_skl, ~, nEval_12] = sklpointsnew_greedy(2, fp, fmin_1, nPart);
[X_skl, ~, nEval_12] = seklpoints_greedyb(2, fp, fmin_1, nPart,initial);
[X_skl, ~, nEval_12] = sejhklpoints_greedy(2, fp, fmin_1, nPart,initial);
[X_skl, ~, nEval_12] = blackboxseklpoints_greedy(2, fp, fmin_1, nPart,initial);

%%%%%%%%%
[X_skl, ~, nEval_12] = blackboxseklpoints_greedyjh(2, fpp, fmin_1, nPart,initial);
%%%%%%%%%
[X_med, ~, nEval_1] = semed_greedyjh(2, fpp, fmin_1, k_1, nPart,initial);
%%%%%%%%%
X_m = gaussmix_rnd(Mu, C, w, nPart);
k = (4 + log(1 + (a - b) * (a - b)')) .^ (-1);
[X_s, ~, nEval] = stein_greedy(2, fp, k, fmin, nPart, []);

t1 = linspace(lb(1), ub(1), nStep)';
t2 = linspace(lb(2), ub(2), nStep)';
T = [repelem(t1, nStep), repmat(t2, nStep, 1)];
p=[];
for i=1:size(T,1)
p(i)= fp(T(i,:));
end
Z = reshape(p, nStep, nStep);




        
        subplot(2,2,1)
contour(t1, t2, Z, 'levelstep', 0.05, 'linewidth', 0.5),title('medpoints');
        
        hold on;
        plot(X_med(:, 1), X_med(:, 2), '.r', 'markersize', 13);
        
subplot(2,2,2)
contour(t1, t2, Z, 'levelstep', 0.05, 'linewidth', 0.5),title('klpoints');
        
        hold on;
        plot(X_kl(:, 1), X_kl(:, 2), '.r', 'markersize', 13);
         text(X_kl(:, 1)+0.02,X_kl(:, 2),num2cell(1:50));
        
      subplot(2,2,3)  
        contour(t1, t2, Z, 'levelstep', 0.05, 'linewidth', 0.5),title('klpointsnew');
        
        hold on;
        plot(X_klnew(:, 1), X_klnew(:, 2), '.r', 'markersize', 13);
        subplot(2,2,4)
        contour(t1, t2, Z, 'levelstep', 0.05, 'linewidth', 0.5),title('medw');
        
        hold on;
        plot(X_wmed(:, 1), X_wmed(:, 2), '.r', 'markersize', 13);
        
        
        
        contour(t1, t2, Z, 'levelstep', 0.05, 'linewidth', 0.5),title('');
        
        hold on;
        plot(X_skl(1:20, 1), X_skl(1:20, 2), '+r', 'markersize', 10);
        hold on;
        for i=21:40
            text(X_skl(i, 1),X_skl(i, 2),num2str(i));
        end
        set(gca,'xtick',0:0.2:1) 
        set(gca,'ytick',0:0.2:1) 
        axis([0 0.2 1,0 0.21])
        plot(X_skl(1:20, 1), X_skl(1:20, 2), '.r', 'markersize', 13);
        hold on;
        plot(X_skl(21:40, 1), X_skl(21:40, 2), '*b', 'markersize', 10);
        
        for i=1:40
            text(X_skl(i, 1),X_skl(i, 2),num2str(i));
        end
        text(X_skl(:, 1)+0.02,X_skl(:, 2),num2cell(1:40));
        
        
        
        contour(t1, t2, Z, 'levelstep', 0.05, 'linewidth', 0.5),title('SMED');
        
        hold on;
        plot(X_med(1:20, 1), X_med(1:20, 2), '.r', 'markersize', 13);
        hold on;
        plot(X_med(21:40, 1), X_med(21:40, 2), '*b', 'markersize', 10);
        text(X_med(:, 1)+0.02,X_med(:, 2),num2cell(1:50));
        
        
        


S=x(1:40, :);
%%%%compute response of these 40 points%%%%
    Y=zeros(40,1);
    for i=1:40
        Y(i)=fp(S(i,:));
    end
    f0=@(x)dace(x,S,Y);
  %%%%%generate the preduce points%%%%%
  xpp=rand(1000,2);
  yxp=zeros(40,1);
  for j=1:40
      yxp(j)=f0(xpp(j,:));
  end
  %%%%%%compute the mse%%%%
  yt=zeros(40,1);
  for j=1:40
      yt(j)=fp(xpp(j,:));
  end
    
  mse = sqrt(sum((yxp - yt) .^ 2)/40);
   0.0226
   0.0255
   0.0223
    0.0163 mse of aklbb40;
  
  
  S=xmm40;
%%%%compute response of these 40 points%%%%
    Y=zeros(40,1);
    for i=1:40
        Y(i)=fp(S(i,:));
    end
    f0=@(x)dace(x,S,Y);
  %%%%%generate the preduce points%%%%%
  xp=rand(1000,5);
  yxp=zeros(40,1);
  for j=1:40
      yxp(j)=f0(xpp(j,:));
  end
  %%%%%%compute the mse%%%%
  yt=zeros(40,1);
  for j=1:40
      yt(j)=fp(xpp(j,:));
  end
    
  mse = sqrt(sum((yxp - yt) .^ 2)/40);
  0.0452
  0.0457
  0.0450
  
  %%%%%smed%%%%%
   S=xp;
%%%%compute response of these 40 points%%%%
    Y=zeros(40,1);
    for i=1:40
        Y(i)=fp(S(i,:));
    end
    f0=@(x)dace(x,S,Y);
  %%%%%generate the preduce points%%%%%
  xp=rand(1000,5);
  yxp=zeros(40,1);
  for j=1:40
      yxp(j)=f0(xpp(j,:));
  end
  %%%%%%compute the mse%%%%
  yt=zeros(40,1);
  for j=1:40
      yt(j)=fp(xpp(j,:));
  end
    
  mse = sqrt(sum((yxp - yt) .^ 2)/40);
  0.0219
  0.0255
   0.0220
   
   
   
   %%%%bijiao jifnezhunze%%%%
   
   S=X_skl;
%%%%compute response of these 40 points%%%%
    Y=zeros(50,1);
    for i=1:50
        Y(i)=fp(S(i,:));
    end
    f0=@(x)dace(x,S,Y);
  %%%%%generate the preduce points%%%%%
  xpp=rand(1000,2);
  yxp=zeros(1000,1);
  for j=1:1000
      yxp(j)=f0(xpp(j,:));
  end
  %%%%%%compute the mse%%%%
  yt=zeros(1000,1);
  for j=1:1000
      yt(j)=fp(xpp(j,:));
  end
    
  e = sum(abs(yxp - yt)/1000);
  e=0.0099 x=bbakl
  mse = sqrt(sum((yxp - yt) .^ 2)/1000);
  mse= 0.0207
  %%%Mmsheji%%%
  S=xmm40;
%%%%compute response of these 40 points%%%%
    Y=zeros(40,1);
    for i=1:40
        Y(i)=fp(S(i,:));
    end
    f0=@(x)dace(x,S,Y);
  %%%%%generate the preduce points%%%%%
  xp=rand(1000,2);
  yxp=zeros(1000,1);
  for j=1:1000
      yxp(j)=f0(xpp(j,:));
  end
  %%%%%%compute the mse%%%%
  yt=zeros(1000,1);
  for j=1:1000
      yt(j)=fp(xpp(j,:));
  end
  
  mse = sqrt(sum((yxp - yt) .^ 2)/1000);
  mse=0.0272
  e = sum(abs(yxp - yt)/1000);
  e=0.0170
  
  
  
  S=xmed;
%%%%compute response of these 40 points%%%%
    Y=zeros(50,1);
    for i=1:50
        Y(i)=fp(S(i,:));
    end
    f0=@(x)dace(x,S,Y);
  %%%%%generate the preduce points%%%%%
  
  
  yxp=zeros(1000,1);
  for j=1:1000
      yxp(j)=f0(xpp(j,:));
  end
  %%%%%%compute the mse%%%%
  yt=zeros(1000,1);
  for j=1:1000
      yt(j)=fp(xpp(j,:));
  end
  mse = sqrt(sum((yxp - yt) .^ 2)/1000);
   0.0033
   
   %%%%%%初开始设计不变的时候
   C=zeros(50,2,50);
   xlcation=[];
yre=[];
for i=1:50
[X_skl, ~, nEval_12] = blackboxseklpoints_greedyjh(2, fpp, fmin_1, nPart,initial);
y=zeros(30,1);
  for j=21:50
      y(j)=fp(X_skl(j,:));
  end
  [m,p]=max(y);
  xlcation(i,:)=X_skl(p,:);
  yre(i)=m;
  C(:,:,i)=X_skl;
end

%%%%%%初始设计是随机的lhd%%%5
Cr=zeros(50,2,50);
   xlcationr=[];
yrer=[];
for i=1:50
    initial=lhd(20,2);
[X_skl, ~, nEval_12] = blackboxseklpoints_greedyjh(2, fpp, fmin_1, nPart,initial);
y=zeros(30,1);
  for j=21:50
      y(j)=fp(X_skl(j,:));
  end
  [m,p]=max(y);
  xlcationr(i,:)=X_skl(p,:);
  yrer(i)=m;
  Cr(:,:,i)=X_skl;
end
%%%%%%%与最优值之间的欧式距离%%%%
xopt=[0.21 0.17];
Xopt=repmat(xopt,50,1);
deklopt=sqrt(sum((Xopt-xlcationr).^2,2));

%%%%%%
 xlcation=[];
yre=[];
for i=1:50

y=zeros(30,1);
  for j=21:50
      y(j)=fp(X_skl(j,:));
  end
  [m,p]=max(y);
  xlcation(i,:)=X_skl(p,:);
  yre(i)=m;
  
end


C=zeros(50,2,100);
 xlcation=[];
yre=[];
for i=1:100
[X_med, ~, nEval_1] = semed_greedyjh(2, fpp, fmin_1, k_1, nPart,initial);
y=zeros(20,1);
  for j=1:20
      y(j)=fp(X_skl(j,:));
  end
  [m,p]=max(y);
  xlcation(i,:)=X_skl(p,:);
  yre(i)=m;
  C(:,:,i)=X_med;
end



%%%%%计算50个kl点预测mse
xpp=rand(1000,2);
mse=zeros(50,1);
for t= 1:50
S=Cr(:,:,t);
%%%%compute response of these 40 points%%%%
    Y=zeros(50,1);
    for i=1:50
        Y(i)=fp(S(i,:));
    end
    f0=@(x)dace(x,S,Y);
  %%%%%generate the preduce points%%%%%
  
  yxp=zeros(1000,1);
  for j=1:1000
      yxp(j)=f0(xpp(j,:));
  end
  %%%%%%compute the mse%%%%
  yt=zeros(1000,1);
  for j=1:1000
      yt(j)=fp(xpp(j,:));
  end
   
  mse(t) = sqrt(sum((yxp - yt) .^ 2)/1000);
end



%%%%%计算kl点预测mse
xpp=rand(1000,2);

S=X_skl;
%%%%compute response of these 40 points%%%%
    Y=zeros(100,1);
    for i=1:100
        Y(i)=fp(S(i,:));
    end
    f0=@(x)dace(x,S,Y);
  %%%%%generate the preduce points%%%%%
  
  yxp=zeros(1000,1);
  for j=1:1000
      yxp(j)=f0(xpp(j,:));
  end
  %%%%%%compute the mse%%%%
  yt=zeros(1000,1);
  for j=1:1000
      yt(j)=fp(xpp(j,:));
  end
   
  mse= sqrt(sum((yxp - yt) .^ 2)/1000);
  mse=0.2475 50 points
  mse=0.0243 100 points;
  %%%%%计算smed的mse
  S=X_med;
    Y=zeros(100,1);
    for i=1:100
        Y(i)=fp(S(i,:));
    end
    f0=@(x)dace(x,S,Y);
  %%%%%generate the preduce points%%%%%
  
  yxp=zeros(1000,1);
  for j=1:1000
      yxp(j)=f0(xpp(j,:));
  end
  %%%%%%compute the mse%%%%
  yt=zeros(1000,1);
  for j=1:1000
      yt(j)=fp(xpp(j,:));
  end
   
  mse= sqrt(sum((yxp - yt) .^ 2)/1000);
  mse=0.3314; 50 points
  mse=0.0188 100 points
  %%%%%计算xmm的mse
  S=xmm;
  S=lhd(100,2);
    Y=zeros(100,1);
    for i=1:100
        Y(i)=fp(S(i,:));
    end
    f0=@(x)dace(x,S,Y);
  %%%%%generate the preduce points%%%%%
  
  yxp=zeros(1000,1);
  for j=1:1000
      yxp(j)=f0(xpp(j,:));
  end
  %%%%%%compute the mse%%%%
  yt=zeros(1000,1);
  for j=1:1000
      yt(j)=fp(xpp(j,:));
  end
   
  mse= sqrt(sum((yxp - yt) .^ 2)/1000);
  mse=0.2737; 
mse=0.4712 lhd 50 points
mse=0.0088 100 lhd
