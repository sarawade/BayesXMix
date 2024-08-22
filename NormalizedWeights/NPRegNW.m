function [w,theta,J,d,cputime,lastState, PP, phi, initialState]=NPRegNW(Y,x,q,p,mcmc,Hyperparameters,varargin)
% function [w,theta,J,d,cputime,lastState]=NPRegNW(Y,x,q,p,MCMC,Hyperparameters,...)
%
% Nonparametric mixture model for regression with normalized weights
% INPUT:
% Y: n x 1 continuous response variable
% x: n x q+p matrix of covariates where the first q columns are discrete
%   and the last p columns are continuous covariates
% q: number of discrete covariates (0 if no discrete covariates)
% p: number of continuous covariates (0 if no continuous covariates)
% mcmc=(burnin, thin, S): 
%   burnin: number of iterations for the burn in period
%   thin: number of iterations for thinning
%   S: MCMC sample size
% Hyperparameters: cell array of hyperparameters for the model
%   (beta0,iC,[alpha1,alpha2]): hyperparameters for the Normal-Inverse Gamma
%       prior for the (beta,sigma)_j
%       beta0 is a q+p+1 x 1 vector
%       iC is a q+p+1 x q+p+1 matrix
%       alpha1,alpha2 are scalars
%   gamma (if q>0): hyperparameters for the Dirichlet prior for the rho_j 
%       gamma is a cell array of length q, where cell h is vector of size
%       Gh, the number of categories of the corresponding covariate.
%   [mu0,c,a1,a2](if p>0): hyperparameters for the Normal-Gamma prior for
%       the (mu_j,tau). They are all p x 1 vectors.
%
% OPTIONAL ARGUMENTS
% Optional arguments must be introduced as 'Property', value. Valid options
% are the following properties:
%   'Prior', the value is a cell array with type of Nonparametric Prior and hyperparameters
%       Type may be:
%           'StickBreaking' (default) for a two parameter stick breaking Prior
%           'Geometric' for Geometric Stick breaking Prior 
%           'DirichletHP' for Dirichlet Process with Gamma hyperprior for
%           the mass parameter
%       Hyperparameters:
%           If Type is 'StickBreaking' or 'Geometric' then the hyperparameters 
%           are (b1,b2). The default is (1,1) corresponding to a Dirichlet 
%           Process Prior with mass parameter 1. 
%           If Type is 'DirichletHP', then the hyperparameters are (b1, b2) 
%           are the parameters for the Gamma hyperprior on the mass parameter 
%           of the Dirichlet process
%   'Phi', the value is the slice sampling constant (scalar). Default value, 0.1
%   'Initial', the value is an initial state of the MCMC sampler
%       (da,ka,Da): initial latent variables
%       rhoa,mua,taua: initial covariate component parameters
%       The default sets ka=0, Da=[], da=(1:n), corresponding to
%       everybody starting in their own cluster with parameters mua=x and all
%       other equal to their prior mean.
%       If the Prior Type is 'DirichletHP', then 'Initial' includes also 
%       massa: initial value for the mass parameter of the Dirichlet Process
%
% OUTPUT:
% w: cell array of length S, contaning sampled posterior weights
% theta=(beta,sigma,rho,mu,tau): posterior samples for the model parameters
% J: S x 1 vector of the number of components sampled at each iteration
% d: S x n matrix of component indices 
% cputime
% lastState: finishing state of the MCMC sampler (to be used as initial
%   state if a larger sample is desired)
%   'PP', the Nonparametric Prior used with hyperparameters.
%       And if Prior Type is 'DirichletHP', PP{3} contains the sequence of
%       sampled mass parameters
%   'phi', the value used as slice sampling constant.
%   'initialState', the  initial state used for the MCMC sampler


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START CPU TIMER %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reversible jump constant
pk=0.5;
% Number of iterations for truncated normal sampling
nitn=10;
% Number of iterations for truncated gamma sampling
nitg=20;
% Sample size
n=length(Y);
% Counter for the number of iterations thinned
ithin= mcmc(2); 
% Counter for saved iterations
it=0;
% Total number of MCMC iterations
MCMC= mcmc(1)+(mcmc(3)-1)*mcmc(2)+1; 
% Augmented X matrix
X=[ones(n,1),x];
% compute C
C=inv(Hyperparameters{2});

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET DEFAULT VALUES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prior for weights
PP=cell(2,1);
PP{1}='StickBreaking';
PP{2}=[1,1];
% Slice sampling constant
phi=0.1;
% Initial State
initialState=cell(3+(q>0)+2*(p>0),1);
% da
initialState{1}=(1:n);
% ka
initialState{2}=zeros(1,n);
% Da
initialState{3}=cell(n,1);
if q==0 % Only continuous variables
    % taua
    initialState{5}=gamrnd(Hyperparameters{4}(:,3),Hyperparameters{4}(:,4).^-1);
    % mua
    initialState{4}=randn(p,n)./kron(sqrt(Hyperparameters{4}(:,2).*initialState{5}),ones(1,n))+kron(Hyperparameters{4}(:,1),ones(1,n));
else % Discrete variables
    % rhoa
    initialState{4}=cell(q,1);
    for h=1:q
        initialState{4}{h}=dirichletrnd(Hyperparameters{4}{h},n);
    end
    clear h
    if p>0 % Continuous variables
        % taua
        initialState{6}=gamrnd(Hyperparameters{5}(:,3),Hyperparameters{5}(:,4).^-1);    
        % mua
        initialState{5}=randn(p,n)./kron(sqrt(Hyperparameters{5}(:,2).*initialState{6}),ones(1,n))+kron(Hyperparameters{5}(:,1),ones(1,n));
    end
end
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OVERRIDE DEFAULT VALUES IF OPTIONAL ARGUMENTS ARE USED %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optarg=length(varargin);
for i=1:optarg/2
    if strcmp(varargin{i*2-1},'Prior')% Prior for weights
        PP=varargin{i*2};
        % Default initial value for mass parameter of Dirichlet process is the mean of the Gamma hyperprior
        if strcmp(PP{1},'DirichletHP')
            initialState{7}=PP{2}(1)/PP{2}(2);
        end
    elseif strcmp(varargin{i*2-1},'Phi')% Slice sampling constant
        phi=varargin{i*2};
    elseif strcmp(varargin{i*2-1},'Initial')% Initial state
        initialState=varargin{i*2};
    else
        disp(['Warning: Invalid property name ',varargin{i*2-1},' default used'])
    end
end
clear i optarg varargin       
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALLOCATION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uniform latent variables allocation
U=cell(n,1);
% Indicators for the uniform latent variables allocation
u=cell(n,1);
% Output variable allocation
J=zeros(mcmc(3),1);
w=cell(mcmc(3),1);
beta=cell(mcmc(3),1);
sigma=cell(mcmc(3),1);
if q>0
    rho=cell(mcmc(3),q);
    if p>0
        mu=cell(mcmc(3),1);
        tau=zeros(mcmc(3),p);
    else
        mu=[];
        tau=[];
    end
else
    rho=[];
    mu=cell(mcmc(3),1);
    tau=zeros(mcmc(3),p);
end
d=zeros(mcmc(3),n);
if strcmp(PP{1},'DirichletHP')
    % Sequence of Dirichlet Process mass parameters when random
    PP{3}=zeros(mcmc(3),1); 
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CURRENT STATE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
da=initialState{1};
ka=initialState{2};
Da=initialState{3};
if q==0 % Only continuous variables
    mua=initialState{4};
    Ja=length(mua);
    taua=initialState{5};
else % Discrete variables
    rhoa=initialState{4};
    Ja=size(rhoa{1},2);
    if p>0 % Continuous variables
        mua=initialState{5};
        taua=initialState{6};
    end
end
if strcmp(PP{1},'DirichletHP')
    massa=initialState{7}; 
end
%clear initialState

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HYPERPARAMETERS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta0=Hyperparameters{1};
iC=Hyperparameters{2};
alpha1=Hyperparameters{3}(1);
alpha2=Hyperparameters{3}(2);
if q==0 % Only continuous variables
    mu0=Hyperparameters{4}(:,1);
    c=Hyperparameters{4}(:,2);
    a1=Hyperparameters{4}(:,3);
    a2=Hyperparameters{4}(:,4);
else % Discrete variables
    gamma=Hyperparameters{4};
    if p>0 % Continuous variables
        mu0=Hyperparameters{5}(:,1);
        c=Hyperparameters{5}(:,2);
        a1=Hyperparameters{5}(:,3);
        a2=Hyperparameters{5}(:,4);
    end
end

%% Start iterations: burn out period %%
for imcmc=1:MCMC
    if mod(imcmc,10)==0
        disp (['Iterations completed= ',int2str(imcmc)]); 
        if strcmp(PP{1},'DirichletHP')
            str=['Mass ', num2str(massa)];
            disp(str)
        end
%        if p>0
%            str = ['Tau ', num2str(taua)'];
%            disp(str)
%        end
        str = ['Max(K)= ', num2str(max(ka)), ', Max(J)= ', num2str(Ja)];
            disp(str)
    end
    %% Update the J (number of components) %%
    if q>0&&p>0
        [Jd,JD,Ja,rhoa,mua]=Jupdate(phi,n,q,p,da,ka,Da,Ja,rhoa,mua,taua,gamma,mu0,c);
    elseif q>0
        [Jd,JD,Ja,rhoa,mua]=Jupdate(phi,n,q,p,da,ka,Da,Ja,rhoa,[],[],gamma,[],[]);
        clear mua;
    else
        [Jd,JD,Ja,rhoa,mua]=Jupdate(phi,n,q,p,da,ka,Da,Ja,[],mua,taua,[],mu0,c);
        clear rhoa;
    end
    
    %% Update the U,u (uniform latent variables and indicators)%%
    if max(ka)>0
        if q>0&&p>0
            [U,u]=Uupdate(x,n,q,p,rhoa,taua,mua,ka,Da);  
        elseif q>0
            [U,u]=Uupdate(x,n,q,p,rhoa,[],[],ka,Da);  
        else
            [U,u]=Uupdate(x,n,q,p,[],taua,mua,ka,Da);  
        end
    end
    
    %% Update Y parameters
    [betaa,sigmaa]=Ypupdate(Y,X,da,Ja,beta0,iC,alpha1,alpha2,C);
    
    %% Update x parameters
    if q>0
        %% Update rho (discrete covariate parameters)
        rhoa=rhoupdate(x(:,1:q),n,q,da,ka,Da,U,u,Ja,rhoa,gamma,nitg);
    end
  
    if p>0
        %% Update mu and tau (continuous covariate parameters)        
        taua=tauupdate(x(:,(q+1):end),n,q,p,da,Ja,mua,taua,ka,Da,U,u,mu0,c,a1,a2,nitg);
        %disp(taua)
        mua=muupdate(x(:,(q+1):end),n,q,p,da,Ja,taua,ka,Da,U,u,mu0,c,nitn);
    end
  
    %% Update w (weights) %%
    if strcmp(PP{1},'DirichletHP')
        massa=massupdate(massa,PP{2},n,da,ka,Da); 
        wa=wupdate(PP{1},massa,Ja,n,da,ka,Da);
    else
        wa=wupdate(PP{1},PP{2},Ja,n,da,ka,Da);
    end
    
    %% Propose label-switching moves %%
    if q>0&&p>0
        [da,Da,sigmaa,betaa,rhoa,mua]=move1(n,q,p,sigmaa,betaa,rhoa,mua,wa,da,Da);
        [da,Da,sigmaa,betaa,rhoa,mua,wa]=move2(n,q,p,sigmaa,betaa,rhoa,mua,wa,da,Da);
    elseif q>0
        [da,Da,sigmaa,betaa,rhoa,mua]=move1(n,q,p,sigmaa,betaa,rhoa,[],wa,da,Da);
        [da,Da,sigmaa,betaa,rhoa,mua,wa]=move2(n,q,p,sigmaa,betaa,rhoa,[],wa,da,Da);
    else
        [da,Da,sigmaa,betaa,rhoa,mua]=move1(n,q,p,sigmaa,betaa,[],mua,wa,da,Da);
        [da,Da,sigmaa,betaa,rhoa,mua,wa]=move2(n,q,p,sigmaa,betaa,[],mua,wa,da,Da);
    end
    
    %% Update d,D (the indices)
    if q>0&&p>0
        [da,Da]=Dupdate(phi,Y,x,n,q,p,sigmaa,betaa,rhoa,mua,taua,wa,Jd,ka,JD);
    elseif q>0
        [da,Da]=Dupdate(phi,Y,x,n,q,p,sigmaa,betaa,rhoa,[],[],wa,Jd,ka,JD);
    else
        [da,Da]=Dupdate(phi,Y,x,n,q,p,sigmaa,betaa,[],mua,taua,wa,Jd,ka,JD);
    end
    
    %% Update k
    if q>0&&p>0
        [ka,Da]=kupdate2(pk,x,n,q,p,Ja,wa,rhoa,taua, mua, ka, Da);
    elseif q>0
        [ka,Da]=kupdate2(pk,x,n,q,p,Ja,wa,rhoa,taua, mua, ka, Da);
    else
        [ka,Da]=kupdate2(pk,x,n,q,p,Ja,wa,[],taua, mua, ka, Da);
    end
    
    %% After burn out period%%

    if imcmc>mcmc(1)
        if ithin==mcmc(2) % Saved
            it=it+1;
            w{it}=wa;
            sigma{it}=sigmaa;
            beta{it}=betaa;
            J(it)=Ja;
            if q>0
                rho{it}=rhoa;
            end
            if p>0
                mu{it}=mua;
                tau(it,:)=taua;
            end
            d(it,:)=da;
            if strcmp(PP{1},'DirichletHP')
                PP{3}(it)=massa; 
            end
            ithin=1;
        else
            ithin=ithin+1;
        end
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE OUTPUT %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Posterior sample of model parameters
theta=cell(2+(q>0)+2*(p>0),1);
theta{1}=beta;
theta{2}=sigma;
if q==0 % Only continuous variables
    theta{3}=mu;
    theta{4}=tau;
else % Discrete variables
    theta{3}=rho;
    if p>0 % Continuous variables
        theta{4}=mu;
        theta{5}=tau;
    end
end
% Last State
lastState=cell(3+(q>0)+2*(p>0),1);
lastState{1}=da;
lastState{2}=ka;
lastState{3}=Da;
if q==0 % Only continuous variables
    lastState{4}=mua;
    lastState{5}=taua;
else % Discrete variables
    lastState{4}=rhoa;
    if p>0 % Continuous variables
        lastState{5}=mua;
        lastState{6}=taua;
    end
end
if strcmp(PP{1},'DirichletHP')
    lastState{7}=massa; 
end
% Time
cputime=toc;
end
