% This article is organized as follows:  
%Section 1 provides the code of generating Kullback-Leibler(KL)points and adaptive Kullback-Leibler (AKL)
% points of a banana-shaped density function described by Haario et al.(1999).
% Section 2 provides the code of generating Kullback-Leibler(KL)points of a
% two-component mixed normal distribution.
% Section 3 provides the code of generating Kullback-Leibler(KL)points of a
% a bivariate normal distribution.
% Section 4 provides the code of generating adaptive Kullback-Leibler (AKL) points and sequential minimum energy design (SMED) points
% of a black box function  described by Fasshauer (2007).

%%%%%% Section 1   %%%%%%%%%
%%%%% The target function: banana shape density function.
fp = @(X)banafuc(X);
d=2; % d is the number of dimensions of the target density.
% Symbolic variables
a = sym('a', [1, 2], 'real');
b = sym('b', [1, 2], 'real');

% Optimisers
lb = [0, 0];
ub = [1, 1];
mu0 = [0, 0];
Sigma0 = diag([25, 25]);
nStep = 100;
fmin_1 = @(f, X)fmin_adares(f, X, d, mu0, Sigma0, 1, 30, lb, ub, [], []);
fmin_2 = @(f, X)fmin_adaresseed(f, X, d, mu0, Sigma0, 1, 30, lb, ub, [], []);

nPart=100;
% Generating 100 MED points of the banana shape density function.
k_1=4*d; % The power parameter of the generalised energy criterion described by Joseph et al (2017).
[X_med, ~, nEval_1] = med_greedy(d, fp, fmin_1, k_1, nPart); 
% Generating 100 KL points of the banana shape density function.
[X_kl, ~, nEval_10] = klpointsbana_greedy(d, fp, fmin_2, nPart);%%


%Generating 50 SMED points and 50 AKL points of the banana shape density
%function, with the initial set is a 20-points MmLHD. 
nPart=50;
%Generating initial set in [0,1]^2£¬which is ontain by an MmLHD include 20
%points.
MmLHD20=importdata('MmLHD20points.txt');
rand('seed',2);
initial=(MmLHD20-rand(20,2))/20;%%%%%%%%%%%


%%%%%%%%Generating the next 30 SMED points based on the initial set
[X_smed, ~, nEval_1] = semed_greedy(2, fp, fmin_2, k_1, nPart,initial);

%%%%%%%%Generating the next 30 AKL points based on the initial set
[X_akl, ~, nEval_12] = aklpoints_greedy(2, fp, fmin_2, nPart,initial);


% Prepare for contour plots.
t1 = linspace(lb(1), ub(1), nStep)';
t2 = linspace(lb(2), ub(2), nStep)';
T = [repelem(t1, nStep), repmat(t2, nStep, 1)];
p=[];
for i=1:size(T,1)
p(i)= fp(T(i,:));
end
Z = reshape(p, nStep, nStep);




        
  %%%%%Scatter graph of KL points for the banana shape function      
contour(t1, t2, Z, 'levelstep', 0.14, 'linewidth', 0.5),title('KL points');  
        hold on;
        plot(X_kl(:, 1), X_kl(:, 2), '.r', 'markersize', 13);
        set(gca,'xtick',[0:0.2:1])
        set(gca,'ytick',[0:0.2:1])
         
        
        
     %%%%%Scatter graph of AKL points and SMED points for the banana shape function   
     
      subplot(1,2,1)
        contour(t1, t2, Z, 'levelstep', 0.14, 'linewidth', 0.5),title('AKL points');
        
        hold on;
        plot(X_akl(1:20, 1), X_akl(1:20, 2), '.r', 'markersize', 13);
        hold on;
        plot(X_akl(21:50, 1), X_akl(21:50, 2), '*b', 'markersize', 9);
        set(gca,'xtick',[0:0.2:1])
        set(gca,'ytick',[0:0.2:1])
        
        subplot(1,2,2)
        contour(t1, t2, Z, 'levelstep', 0.14, 'linewidth', 0.5),title('SMED points');
        
        hold on;
        plot(X_smed(1:20, 1), X_smed(1:20, 2), '.r', 'markersize', 13);
        hold on;
        plot(X_smed(21:50, 1), X_smed(21:50, 2), '*b', 'markersize', 9);
        set(gca,'xtick',[0:0.2:1])
       set(gca,'ytick',[0:0.2:1])
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Section 2 %%%%%%%%%%%%%%%%%%%%%%%%
% The target is a two-component mixed normal distribution with the mean
% vector mu_1=(-1.5,0)',mu_2=(1.5,0)' and covariance matrix
% Sigma_1=Sigma_2=I,where I is a identity matrix.

% Target density  functions
[Mu, C, w] = gmparam_2c();
fp = @(X)fp_gaussmix(X, Mu, C, w);

% Symbolic variables
a = sym('a', [1, 2], 'real');
b = sym('b', [1, 2], 'real');
d=2; % the dimension of target density function.

% Optimisers
lb = [-5, -5];
ub = [5, 5];
mu0 = [0, 0];
Sigma0 = diag([16, 16]);
nStep = 100;
fmin_2 = @(f, X)fmin_adaresseed(f, X, d, mu0, Sigma0, 1, 30, lb, ub, [], []);


% Extensible point sets
nPart = 100;


% (1) Monte Carlo point set
X_m = gaussmix_rnd(Mu, C, w, nPart);

% (2) MED point set
nEval = ones(nPart, 1);
X_1 = cumset(X_m);
nEval_1 = log(cumsum(nEval));
k_1=4*d;
[X_med, ~, nEval_1] = med_greedy(d, fp, fmin_1, k_1, nPart);

% (3) KL point set
[X_kl, ~, nEval_2] = klpointsmixnr_greedy(d, fp, fmin_2, nPart);



% Prepare for contour plots.
t1 = linspace(lb(1), ub(1), nStep)';
t2 = linspace(lb(2), ub(2), nStep)';
T = [repelem(t1, nStep), repmat(t2, nStep, 1)];
p = fp(T);
Z = reshape(p, nStep, nStep);


 %%%%% Scatter graph of MC points, MED points and KL points for the
 %%%%% mixed normal distribution.
contour(t1, t2, Z, 'levelstep', 0.007, 'linewidth', 0.5),title('MC points');
        
        hold on;
        plot(X_m(:, 1), X_m(:, 2), '.r', 'markersize', 10);
        set(gca,'xtick',[-5:2:5])
        set(gca,'ytick',[-5:2:5])
        

contour(t1, t2, Z, 'levelstep', 0.007, 'linewidth', 0.5),title('MED points');
        
        hold on;
        plot(X_med(:, 1), X_med(:, 2), '.r', 'markersize', 10);
        set(gca,'xtick',[-5:2:5])
        set(gca,'ytick',[-5:2:5])

        
contour(t1, t2, Z, 'levelstep', 0.007, 'linewidth', 0.5),title('KL points');
        
        hold on;
        plot(X_kl(:, 1), X_kl(:, 2), '.r', 'markersize', 10);
        set(gca,'xtick',[-5:2:5])
        set(gca,'ytick',[-5:2:5])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
%%%%%%%%%%%%% Section 3 %%%%%%%%%%%%%%%%%%%%%%%%    
% The target is a bivariate normal distribution.

% Target density. 
% (1) The target is a bivariate normal distribution N(0,I), where I is a 2*2
% identity matrix.
[Mu, C, w] = gmparam_2c;
Mu=[0 0;0 0];
w=[1;0];
fp = @(X)fp_gaussmix(X, Mu, C, w);


% (2) The target is a bivariate normal distribution N(0,R^2), where R is a 2*2
% matrix with the R(i,j)=0.5^{|i-j|};
[Mu, C, w] = gmparam_2c
Mu=[0 0;0 0];
C(:,:,1)=[1 0.5;0.5 1];
C(:,:,2)=[1 0.5;0.5 1];
w=[1;0];
fp = @(X)fp_gaussmix(X, Mu, C, w);

% Symbolic variables
d=2; % the number of dimensions of the target density.
a = sym('a', [1, d], 'real');
b = sym('b', [1, d], 'real');


% Optimisers
lb = [-5, -5];
ub = [5, 5];
mu0 = [0, 0];
Sigma0 = diag([16, 16]);
nStep = 100;
fmin_1 = @(f, X)fmin_adares(f, X, d, mu0, Sigma0, 1, 30, lb, ub, [], []);


% Extensible point sets for the target distribution N(0,I)
nPart = 100;
% (1) Monte Carlo point aet
X_m = gaussmix_rnd(Mu, C, w, nPart);

% (2) MED point set 
k_1=4*d;
[X_med, ~, ~] = med_greedy(d, fp, fmin_1, k_1, nPart);

% (3) KL point set 
[X_kl, ~, ~] = klpointsnriid_greedy(d, fp, fmin_1, nPart);


% Extensible point sets for the target distributionN(0,R^2), where R is a 2*2
% matrix with the R(i,j)=0.5^{|i-j|};
nPart = 100;
% (1) Monte Carlo point aet
X_m = gaussmix_rnd(Mu, C, w, nPart);

% (2) MED point set 
k_1=4*d;
[X_med, ~, ~] = med_greedy(d, fp, fmin_1, k_1, nPart);

% (3) KL point set 
[X_kl, ~, ~] = klpointsnr05_greedy(d, fp, fmin_1, nPart);



% Prepare for contour plots
t1 = linspace(lb(1), ub(1), nStep)';
t2 = linspace(lb(2), ub(2), nStep)';
T = [repelem(t1, nStep), repmat(t2, nStep, 1)];
p=[];
for i=1:size(T,1)
p(i)= fp(T(i,:));
end
Z = reshape(p, nStep, nStep);



 %%%%%Scatter graph of MC points, MED points and KL points for the
 %%%%%bivariate normal distribution.
contour(t1, t2, Z, 'levelstep', 0.014, 'linewidth', 0.5),title('MC points');
        
        hold on;
        plot(X_m(:, 1), X_m(:, 2), '.r', 'markersize', 10);
        set(gca,'xtick',-5:2:5)
        set(gca,'ytick',-5:2:5)
        
        
contour(t1, t2, Z, 'levelstep', 0.014, 'linewidth', 0.5),title('MED points');
        
        hold on;
        plot(X_med(:, 1), X_med(:, 2), '.r', 'markersize', 10);
        set(gca,'xtick',-5:2:5)
        set(gca,'ytick',-5:2:5)
        
        
        
        
contour(t1, t2, Z, 'levelstep', 0.014, 'linewidth', 0.5),title('KL points');
        
        hold on;
        plot(X_kl(:, 1), X_kl(:, 2), '.r', 'markersize', 10);
        set(gca,'xtick',-5:2:5)
        set(gca,'ytick',-5:2:5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%%%%%%%%%%%%% Section 4 %%%%%%%%%%%%%%%%%%%%%%%%  
%%%The target black box function.
fp = @(X)bb(X);

% Symbolic variables
d=2;% The number of dimensions of the target density.
a = sym('a', [1, 2], 'real');
b = sym('b', [1, 2], 'real');

% Optimisers
lb = [0, 0];
ub = [1, 1];
mu0 = [0, 0];
Sigma0 = diag([25, 25]);
nStep = 100;

fmin_1 = @(f, X)fmin_adares(f, X, d, mu0, Sigma0, 1, 30, lb, ub, [], []);
fmin_4 = @(f, X)fmin_adaresseed8(f, X, d, mu0, Sigma0, 1, 30, lb, ub, [], []);


% Extensible point sets
nPart = 40;
% Generating initial set in [0,1]^2£¬which is obntain by an MmLHD include 20
% points.
MmLHD20=importdata('MmLHD20points.txt');
rand('seed',2);
initial=(MmLHD20-rand(20,2))/20;


%%%%%%%%%%%%% Generating 20 SMED points with the initial point set
%%%%%%%%%%%%% described by Joseph et al (2017).
k_1=4*d;% k is the exponent of the energy criterion which used to generate sequential minimum energy design (SMED).
[X_smed, ~, nEval_1] = bbsemed_greedy(d, fp, fmin_4, k_1, nPart,initial);

%%%%%%%%%%%%% Generating the next 20 AKL points based on the same initial set.
[X_akl, ~, nEval_12] = blackboxseklpoints_greedy(d, fp, fmin_4, nPart,initial);


% Prepare for contour plots.
t1 = linspace(lb(1), ub(1), nStep)';
t2 = linspace(lb(2), ub(2), nStep)';
T = [repelem(t1, nStep), repmat(t2, nStep, 1)];
p=[];
for i=1:size(T,1)
p(i)= fp(T(i,:));
end
Z = reshape(p, nStep, nStep);

  
 %%%%%Scatter graph of AKL points and SMED points for the black box
 %%%%%function.
 
 x=[0.21,0.17]; % x is the global optimum of the target black box function.
 
        subplot(1,2,1)
         [c,h]=contourf(t1, t2, Z, 11, 'levelstep', 0.05, 'linewidth', 0.5),title('AKL points');
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
plot(X_akl(1:20, 1), X_akl(1:20, 2), '+r', 'markersize', 11);
hold on;
for i=21:40
text(X_akl(i, 1),X_akl(i, 2),num2str(i),'Color','blue');
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
plot(X_smed(1:20, 1), X_smed(1:20, 2), '+r', 'markersize', 11);
hold on;
for i=21:40
text(X_smed(i, 1),X_smed(i, 2),num2str(i),'Color','blue');
end
set(gca,'xtick',0:0.2:1)
set(gca,'ytick',0:0.2:1)
 hold on;
plot(x(1, 1), x(1, 2), '*m', 'markersize', 10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
