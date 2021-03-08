[Mu, C, w] = gmparam_2c();
fp = @(X)frifunc(X);
fu = @(X)fscr_gaussmix(X, Mu, C, w);
fp = @(X)expfrifunc(X);
% Symbolic variables
a = sym('a', [1, 5], 'real');
b = sym('b', [1, 5], 'real');

% Optimisers
lb = [0, 0,0,0,0];
ub = [1, 1,1,1,1];
mu0 = [0, 0,0,0,0];
Sigma0 = diag([25, 25,25,25,25]);
nStep = 100;
fmin_1 = @(f, X)fmin_adares(f, X, 5, mu0, Sigma0, 1, 30, lb, ub, [], []);
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
[X_kl, ~, nEval_10] = klpoints_greedy(2, fp, fmin_1, nPart);
[X_kl, ~, nEval_10] = klpointsk_greedy(2, fp, fmin_1, nPart);
[X_klnew, ~, nEval_11] = klpointsnew_greedy(2, fp, fmin_1, nPart);
[X_skl, ~, nEval_12] = sklpoints_greedy(2, fp, fmin_1, nPart);
[X_skl, ~, nEval_12] = sklpointsnew_greedy(2, fp, fmin_1, nPart);
[X_skl, ~, nEval_12] = seklpoints_greedyb(2, fp, fmin_1, nPart,initial);
[X_skl, ~, nEval_12] = sejhklpoints_greedy(2, fp, fmin_1, nPart,initial);
[X_skl, ~, nEval_12] = blackboxseklpoints_greedy(2, fp, fmin_1, nPart,initial);
%%%%%%%
[X_skl, ~, nEval_12] = seklpoints_greedylcvxx(5, fp, fmin_1, nPart,initial);
%%%%%%

[X_skl, ~, nEval_12] = seklpoints_greedylcvxd(5, fp, fmin_1, nPart,initial);
X_m = gaussmix_rnd(Mu, C, w, nPart);
k = (4 + log(1 + (a - b) * (a - b)')) .^ (-1);
[X_s, ~, nEval] = stein_greedy(2, fp, k, fmin, nPart, []);

t1 = linspace(lb(1), ub(1), nStep)';
t2 = linspace(lb(2), ub(2), nStep)';
t3 = linspace(lb(2), ub(2), nStep)';
t4 = linspace(lb(2), ub(2), nStep)';
t5 = linspace(lb(2), ub(2), nStep)';
T = [repelem(t1, nStep), repmat(t2, nStep, 1),repelem(t3, nStep),repelem(t4, nStep),repelem(t5, nStep)];

T = [repelem(t1, nStep), repmat(t2, nStep, 1),repelem(t3, nStep)];
p=[];
for i=1:size(T,1)
p(i)= fp(T(i,:));
end
Z = reshape(p, nStep, nStep);




        
        subplot(2,2,1)
contour(t1, t2, Z, 'levelstep', 0.5, 'linewidth', 0.5),title('medpoints');
        
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
        
        
        
        contour(t1, t2, Z, 'levelstep', 0.7, 'linewidth', 0.5),title('');
        
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
        
        
        
        contour(t3, t4, Z, 'levelstep', 0.7, 'linewidth', 0.5),title('');
        
        hold on;
        plot(X_skl(1:20, 1), X_skl(1:20, 2), '+r', 'markersize', 10);
        hold on;
        for i=21:40
            text(X_skl(i, 1),X_skl(i, 2),num2str(i));
        end
        
        contour(t1, t2, Z, 'levelstep', 0.7, 'linewidth', 0.5),title('');
        
        hold on;
        plot(xkl01(1:20, 1), xkl01(1:20, 2), '+r', 'markersize', 10);
        hold on;
        for i=21:40
            text(xkl01(i, 1),xkl01(i, 2),num2str(i));
        end
        
        contour(t1, t2, Z, 'levelstep', 0.05, 'linewidth', 0.5),title('SMED');
        
        hold on;
        plot(X_med(1:20, 1), X_med(1:20, 2), '.r', 'markersize', 13);
        hold on;
        plot(X_med(21:40, 1), X_med(21:40, 2), '*b', 'markersize', 10);
        text(X_med(:, 1)+0.02,X_med(:, 2),num2cell(1:50));
        
        
        subplot(1,2,1)
         [c,h]=contourf(t1, t2, Z, 11, 'levelstep', 0.05, 'linewidth', 0.5),title('SKL points');
map = [0.52 0.8 0.92
0.8 0.7 0.8
0.98 0.94 0.90
1 0.75 0.79
0.94 1 0.94];
colormap(map);
set(h,'ShowText','on','LevelList',[-0.1 0 0.1 0.2 0.3 0 .4 0.5 0.6 0.7 0.8 0.9]);
clabel(c,h,'LabelSpacing',10000);
set(gca,'xtick',0:0.2:1)
set(gca,'ytick',0:0.2:1)
hold on;
plot(X_skl(1:20, 1), X_skl(1:20, 2), '+r', 'markersize', 11);
hold on;
for i=21:40
text(X_skl(i, 1),X_skl(i, 2),num2str(i),'Color','blue');
end
set(gca,'xtick',0:0.2:1)
set(gca,'ytick',0:0.2:1)
 hold on;
plot(x(1, 1), x(1, 2), '*m', 'markersize', 10);
        
     subplot(1,2,2)   
        [c,h]=contourf(t1, t2, Z, 11, 'levelstep', 0.05, 'linewidth', 0.5),title('SMED points');
map = [0.52 0.8 0.92
0.8 0.7 0.8
0.98 0.94 0.90
1 0.75 0.79
0.94 1 0.94];
colormap(map);
set(h,'ShowText','on','LevelList',[-0.1 0 0.1 0.2 0.3 0 .4 0.5 0.6 0.7 0.8 0.9]);
clabel(c,h,'LabelSpacing',10000);
set(gca,'xtick',0:0.2:1)
set(gca,'ytick',0:0.2:1)
hold on;
plot(X_med(1:20, 1), X_med(1:20, 2), '+r', 'markersize', 11);
hold on;
for i=21:40
text(X_med(i, 1),X_med(i, 2),num2str(i),'Color','blue');
end
set(gca,'xtick',0:0.2:1)
set(gca,'ytick',0:0.2:1)
 hold on;
plot(x(1, 1), x(1, 2), '*m', 'markersize', 10);
      
        
        set(gcf, 'renderer', 'painters');
set(gcf, 'units', 'centimeters');
set(gcf, 'position', [3, 3, 10, 17.5]);

% Print setting
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'paperpositionmode', 'manual');
set(gcf, 'paperposition', [1, 4, 30, 23.5]);
set(gcf, 'papertype', '<custom>');
set(gcf, 'papersize', [29, 40.5]);

% Print to PDF
name = sprintf('banafunc4meth');
print(name, '-dpdf');


       [c,h]=contourf(t1, t2, Z, 11, 'levelstep', 0.05, 'linewidth', 0.5),title('');
map = [0.52 0.8 0.92
0.8 0.7 0.8
0.98 0.94 0.90
1 0.75 0.79
0.94 1 0.94];
colormap(map);
set(h,'ShowText','on','LevelList',[-0.1 0 0.1 0.2 0.3 0 .4 0.5 0.6 0.7 0.8 0.9]);
clabel(c,h,'LabelSpacing',10000);
set(gca,'xtick',0:0.2:1)
set(gca,'ytick',0:0.2:1)
hold on;
plot(x(1, 1), x(1, 2), '*r', 'markersize', 10);



 [c,h]=contourf(t1, t2, Z, 11, 'levelstep', 0.05, 'linewidth', 0.5),title('');
map = [0.52 0.8 0.92
0.8 0.7 0.8
0.98 0.94 0.90
1 0.75 0.79
0.94 1 0.94];
colormap(map);
set(h,'ShowText','on','LevelList',[-0.1 0 0.1 0.2 0.3 0 .4 0.5 0.6 0.7 0.8 0.9]);
clabel(c,h,'LabelSpacing',10000);
set(gca,'xtick',0:0.2:1)
set(gca,'ytick',0:0.2:1)
hold on;
plot(X_med(:, 1), X_med(:, 2), '*r', 'markersize', 10);


%%%%%%%%用01区间上的40个kl点作预测%%%%%
S=xkl0189(1:40, :);
%%%%compute response of these 40 points%%%%
    Y=zeros(40,1);
    for i=1:40
        Y(i)=fp(S(i,:));
    end
    f0=@(x)dace(x,S,Y);
  %%%%%generate the preduce points%%%%%
  xp=rand(1000,5);
  yxp=zeros(1000,1);
  for j=1:1000
      yxp(j)=f0(xp(j,:));
  end
  %%%%%%compute the mse%%%%
  yt=zeros(1000,1);
  for j=1:1000
      yt(j)=fp(xp(j,:));
  end
    
  mse = sqrt(sum((yxp - yt) .^ 2)/1000);
  mse=3.3394;
  
  %%%%%%%%40ge maximin points
  S=xmm(1:40, :);
%%%%compute response of these 40 points%%%%
    Y=zeros(40,1);
    for i=1:40
        Y(i)=fp(S(i,:));
    end
    f0=@(x)dace(x,S,Y);
  %%%%%generate the preduce points%%%%%
  xp=rand(1000,5);
  yxp=zeros(1000,1);
  for j=1:1000
      yxp(j)=f0(xp(j,:));
  end
  %%%%%%compute the mse%%%%
  yt=zeros(1000,1);
  for j=1:1000
      yt(j)=fp(xp(j,:));
  end
    
  mse = sqrt(sum((yxp - yt) .^ 2)/1000);
  mse=2.3226;
  
  %%%%%%%%用01区间上的后20个kl点作预测%%%%%
S=xkl0189(21:40, :);
%%%%compute response of these 40 points%%%%
    Y=zeros(20,1);
    for i=1:20
        Y(i)=fp(S(i,:));
    end
    f0=@(x)dace(x,S,Y);
  %%%%%generate the preduce points%%%%%
  xp=rand(1000,5);
  yxp=zeros(1000,1);
  for j=1:1000
      yxp(j)=f0(xp(j,:));
  end
  %%%%%%compute the mse%%%%
  yt=zeros(1000,1);
  for j=1:1000
      yt(j)=fp(xp(j,:));
  end
    
  mse = sqrt(sum((yxp - yt) .^ 2)/1000);
  mse= 7.1401;
  
   %%%%%%%%用01区间上的前20个kl点作预测%%%%%
S=xkl01in(1:20, :);
%%%%compute response of these 40 points%%%%
    Y=zeros(20,1);
    for i=1:20
        Y(i)=fp(S(i,:));
    end
    f0=@(x)dace(x,S,Y);
  %%%%%generate the preduce points%%%%%
  xp=rand(1000,5);
  yxp=zeros(1000,1);
  for j=1:1000
      yxp(j)=f0(xp(j,:));
  end
  %%%%%%compute the mse%%%%
  yt=zeros(1000,1);
  for j=1:1000
      yt(j)=fp(xp(j,:));
  end
    
  mse = sqrt(sum((yxp - yt) .^ 2)/1000);
  mse=3.4626;
  
  
  S=xkl01(21:40, :);
%%%%compute response of these 40 points%%%%
    Y=zeros(20,1);
    for i=1:20
        Y(i)=fp(S(i,:));
    end
    
    S=xkl01(1:20, :);
%%%%compute response of these 40 points%%%%
    YI=zeros(20,1);
    for i=1:20
        YI(i)=fp(S(i,:));
    end