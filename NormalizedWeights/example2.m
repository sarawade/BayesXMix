%% EXAMPLE %%
% Example 2 from Wade and Inacio (2024): demonstrating the drawbacks of the
% conditional approach with single weights

% Clear workspace %
% Comment this line if clearing is not wanted
clear

% Set random number generator seed for data simulation %
% Comment this lines if seed should not be fixed
seed=32455;
rng(seed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=csvread('ex2/ex2train.csv',1,1);
Y=M(:,2);
x=M(:,1);
m_true = M(:,3);
clear M

% Sample size
n=size(Y,1);
% Number of discrete covariates
q=0;
% Number of continuous covariates
p=1;

% Define new covariate values
M=csvread('ex2/ex2test.csv',1,1);
x_new = M(:,1);
n_new=size(x_new,1);
% For simulated data, we can calculate the true predictive mean for
% comparison with estimates
m_pred_true = M(:,2);
clear M

% Define Y points at which predictive density will be estimated
Y_grid=(-1:.02:10)';
ny=length(Y_grid);
% For simulated data, we can calculate the true predictive density for
% comparison with estimates
Y_fpred_true=zeros(n_new,ny);
Y_fpred_true(x_new<=2,:)=normpdf(repmat(Y_grid',sum(x_new<=2),1),repmat(m_pred_true(x_new<=2),1,ny),0.2);
Y_fpred_true(x_new>2 & x_new<=5,:)=normpdf(repmat(Y_grid',sum(x_new>2 & x_new<=5),1),repmat(m_pred_true(x_new>2 & x_new<=5),1,ny),0.05);
Y_fpred_true(x_new>5,:)=normpdf(repmat(Y_grid',sum(x_new>5),1),repmat(m_pred_true(x_new>5),1,ny),repmat((x_new(x_new>5)-5).^2/15+0.01,1,ny));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE PARAMETERS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MCMC PARAMETERS
% MCMC Sample size
S=5000;
% Burn in period
burnin=5000;
% Thinning
thin=1;

% Generate data structure for NPRegNW function and delete auxiliary
% variables
mcmc=[burnin,thin,S];
clear S burnin thin

%% Base measure hyperparameters

% Hyperparameters for the Normal-InverseGamma prior for the (beta,sigma^2)_j
% beta0 is the mean for beta|sigma^2
beta0=[0.3256765, 0.7272275]';
% sigma^2*iC is the variance for beta|sigma^2
iC=diag([100,10]);
% alpha1/alpha2 is the mean for 1/sigma^2
alpha1=2;
alpha2=1/10;

% Hyperparameters for the Normal-Gamma prior for the (mu_j,tau)
% mu0 is the mean for mu|tau
mu0=[3.927097]';
% c*tau is the precision for mu|tau 
c=[1/10]';
% a1/a2 is the mean for tau
a1=[2]';
a2=[1]';

% Generate data structure for NPRegNW function and clear auxiliary
% variables
Hyperparameters=cell(5,1);
Hyperparameters{1}=beta0;
Hyperparameters{2}=iC;
Hyperparameters{3}=[alpha1,alpha2];
%Hyperparameters{4}=gamma;
Hyperparameters{4}=[mu0,c,a1,a2];

%% OPTIONAL PARAMETERS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If undefined, default values will be used

%% Slice sampling constant
% Default value is 0.1
% phi=0.01;

%% Nonparametric prior specification
% First argument is the type. May be 'StickBreaking' or 'Geometric'
% Second argument is Hyperparameters.
% The default is 'StickBreaking' with [1,1]
PP=cell(2,1);
PP{1}='StickBreaking';
PP{2}=[1,1];
% Example 2
% PP{1}='Geometric';
% PP{2}=[1,2];
% Example 3
% PP{1}='DirichletHP';
% PP{2}=[.5,.5];
% massa=1;

%% Initial state: randomized alternative
% Default starts each with one mixture component for each observation,
% with covariate parameters sampled from the prior and empty latent
% variables

% Indices
cl=3; % Number of random initial clusters
da=discreternd((1:cl)',n)'; % only one component with every observation associated to it

% Latent model variables
ka=5*ones(n,1); % No latent variables
Da=cell(n,1); 
for i=1:n
    Da{i}=discreternd((1:2*cl)',ka(i))';
end

% Parameters for continuous covariates
taua=a1./a2; % precision set at prior mean
mua=randn(1,2*cl)./sqrt(taua.*c)+mu0;


% Generate data structure for NPRegNW function and clear auxiliary
% variables
initialState=cell(5,1);
initialState{1}=da;
initialState{2}=ka;
initialState{3}=Da;
initialState{4}=mua;
initialState{5}=taua;
%if strcmp(PP{1},'DirichletHP')
%    initialState{6}=massa;
%    clear massa
%end
clear h da ka Da rhoa mua taua G h beta0 iC alpha1 alpha2 gamma mu0 c a1 a2

%% POSTERIOR SAMPLING %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The first six imput elements are required.
%[w,theta,J,d,cputime,lastState]=NPRegNW(Y,x,q,p,mcmc,Hyperparameters);

% The rest are optional, with format 'Property name', Property value
% [w,theta,J,d,cputime,lastState]=NPRegNW(Y,x,q,p,mcmc,Hyperparameters,'Initial',initialState, 'Prior',PP,'Phi',phi);
% The additional values can be returned by the algorithm to avoid confusion
% about the values used
%[w,theta,J,d,cputime,lastState,PP,phi,initialState]=NPRegNW(Y,x,q,p,mcmc,Hyperparameters,'Initial',initialState, 'Prior',PP,'Phi',phi);

%Example
[w,theta,J,d,cputime,lastState,PP,phi,initialState]=NPRegNW(Y,x,q,p,mcmc,Hyperparameters,'Initial',initialState);

% To continue MCMC iterations, set initialState=lastState
% mcmc=[3000,1,5000];
% [w,theta,J,d,cputime,lastState,PP]=NPRegNW(Y,x,q,p,mcmc,Hyperparameters,'Prior',PP,'Initial',lastState,'Phi',phi);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST-PROCESSING %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inference for configurations
[CS,dR,pdR,last_s_dR]=relabel(d);
S=mcmc(3);

%% Inference for prediction (mean functional)
 
% Probability for pointwise credible intervals
% Optional argument, use only if intervals are required 
cred_prob=0.95;

[Y_pred,Y_pred_mat,Y_pred_CI]=NPRegNWprediction(x_new,q,p,w,theta,J,cred_prob);
% If no credible intervals are needed, use
% [Y_pred]=NPRegNWprediction(x_new,q,p,w,theta,J);

% Compute L2 prediction error
l2_err_pred=sqrt(sum((m_pred_true-Y_pred).^2)/n_new);

ec_pred = sum((m_pred_true>=Y_pred_CI(:,1)) & (m_pred_true<=Y_pred_CI(:,2)))/n_new;

ci_length = mean(Y_pred_CI(:,2)-Y_pred_CI(:,1));
    
%% Predictive density estimates

[Y_fpred,l_fpred,u_fpred]=NPRegNWfprediction(Y_grid,x_new,q,p,w,theta,J,0.95);

% Compute L_1 error between true and estimated (mean posterior) 
% conditional density for each new covariate value 
L1_f_pred=sum(abs(Y_fpred_true-Y_fpred),2)*(Y_grid(2)-Y_grid(1));
% Simple average
savgL1_f_pred=mean(L1_f_pred);

%% SAVE WORKSPACE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save 'ex2/example2_nwreg_v3'

csvwrite('ex2/ex2_pred_nwreg_v3.csv',[Y_pred,Y_pred_CI(:,1:2)]);
csvwrite('ex2/ex2_fpred_nwreg_v3.csv',Y_fpred);
csvwrite('ex2/ex2_lfpred_nwreg_v3.csv',l_fpred);
csvwrite('ex2/ex2_ufpred_nwreg_v3.csv',u_fpred);
csvwrite('ex2/ex2_config_nwreg_v3.csv',d);

