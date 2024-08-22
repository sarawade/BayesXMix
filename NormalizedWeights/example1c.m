%% EXAMPLE %%
% Demonstrating how to use the NPRegNW.m function for implementing the
% Nonparametric regression model with normalized weights.

% Clear workspace %
% Comment this line if clearing is not wanted
clear

% Set random number generator seed for data simulation %
% Comment this lines if seed should not be fixed
%seed=32455;
%rng(seed);
%s = RandStream('mcg16807','Seed',seed);
%RandStream.setDefaultStream(s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=csvread('ex1train_mmerrors.csv',1,1);
Y=M(:,3);
x=M(:,1:2);
% True regression function
m_true = 5-log(x(:,1)+2);
clear M

% Sample size
n=size(Y,1);
% Number of discrete covariates
q=0;
% Number of continuous covariates
p=2;


% Define new covariate values
x_new=csvread('ex1test_mmerrors.csv',1,1);
n_new=size(x_new,1);
% For simulated data, we can calculate the true predictive mean for
% comparison with estimates
m_pred_true=5-log(x_new(:,1)+2);

% Define Y points at which predictive density will be estimated
Y_grid=(2.3:.02:5.4)';
ny=length(Y_grid);
% For simulated data, we can calculate the true predictive density for
% comparison with estimates
%Y_fpred_true=normpdf(repmat(Y_grid',n_new,1),repmat(m_pred_true,1,ny),0.05);
Y_fpred_true = csvread('ex1test_mmerrors_dens.csv',1,1);


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
beta0=[4.1064, -0.21736, 0.07534]'; 
% sigma^2*iC is the variance for beta|sigma^2
iC=diag([100,4,4]);
% alpha1/alpha2 is the mean for 1/sigma^2
alpha1=2;
alpha2=1/100;

% Hyperparameters for the Normal-Gamma prior for the (mu_j,tau)
% mu0 is the mean for mu|tau
mu0=[3.477661, 1.086]';
% c*tau is the precision for mu|tau 
c=[1/40,1/40]';
% a1/a2 is the mean for tau
a1=[2,2]';
a2=[0.1,0.1]';

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

%% Initial state
% Default starts each with one mixture component for each observation,
% with covariate parameters sampled from the prior and empty latent
% variables

% Indices
%da=(1:n); % n components with observation i associated to component i
da=ones(1,n); % only one component with every observation associated to it

% Latent model variables
ka=zeros(n,1); % No latent variables
Da=cell(n,1); % Empty array (because ka_i is zero for all i) of adequate size

% Parameters for discrete covariates
rhoa=cell(q,1);
% Initalized at the prior mean
for h=1:q
    %rhoa{h}=kron(ones(1,n),gamma{h}/sum(gamma{h})); % n components
    rhoa{h}=gamma{h}/sum(gamma{h}); % 1 component
end

% Parameters for continuous covariates
%mua=x(:,q+1:end)'; % n components with means set at data points
mua=sum(x(:,q+1:end))'/n; % one component with mean set at sample mean
taua=a1./a2; % precision set at prior mean

% Generate data structure for NPRegNW function and clear auxiliary
% variables
initialState=cell(6,1);
initialState{1}=da;
initialState{2}=ka;
initialState{3}=Da;
initialState{4}=rhoa;
initialState{5}=mua;
initialState{6}=taua;
if strcmp(PP{1},'DirichletHP')
    initialState{7}=massa;
    clear massa
end
clear h da ka Da rhoa mua taua G h beta0 iC alpha1 alpha2 gamma mu0 c a1 a2

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
% [w,theta,J,d,cputime,lastState,PP,phi,initialState]=NPRegNW(Y,x,q,p,mcmc,Hyperparameters,'Initial',initialState, 'Prior',PP);

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
% % Define new covariates
% x_new=[(-6:.1:6)';(-6:.1:6)'];
% n_new=size(x_new,1);
% x_new=[[zeros(n_new/2,1);ones(n_new/2,1)],x_new];
% Probability for pointwise credible intervals
% Optional argument, use only if intervals are required 
cred_prob=0.95;

[Y_pred,Y_pred_mat,Y_pred_CI]=NPRegNWprediction(x_new,q,p,w,theta,J,cred_prob);
% If no credible intervals are needed, use
% [Y_pred]=NPRegNWprediction(x_new,q,p,w,theta,J);

%For the data
%[Y_mean,Y_mean_mat,Y_mean_CI]=NPRegNWprediction(x,q,p,w,theta,J,cred_prob);
%Y_mean_true=5./(1+exp(-x));

% Mean square error
%MSE=sum(((Y_mean-Y).^2))/n;

% Compute L2 prediction error
% l2_err=sqrt(sum((Y_pred_true-(sum(Y_pred,1)/mcmc(1))').^2)/n_new);
l2_err_pred=sqrt(sum((m_pred_true-Y_pred).^2)/n_new);
%l2_err_pred_median=sqrt(sum((Y_pred_true-Y_pred_CI(:,3)).^2)/n_new);
% Compute L1 prediction error
% l1_err=sum(abs(Y_pred_true-(sum(Y_pred,1)/mcmc(1))'))/n_new;
%l1_err_pred=sum(abs(Y_pred_true-Y_pred))/n_new;
%l1_err_pred_median=sum(abs(Y_pred_true-Y_pred_CI(:,3)))/n_new;
% For the data
%l2_err_mean=sqrt(sum((Y_mean_true-Y_mean).^2)/n);
%l2_err_mean_median=sqrt(sum((Y_mean_true-Y_mean_CI(:,3)).^2)/n);
%l1_err_mean=sum(abs(Y_mean_true-Y_mean))/n;
%l1_err_mean_median=sum(abs(Y_mean_true-Y_mean_CI(:,3)))/n;


% Smoothed histograms representing estimated posterior density for the
% regression mean
%Y_pred_hist=zeros(100,n_new);
%Y_mat=zeros(100,n_new);
%for i=1:n_new
%    [Y_pred_hist(:,i),Y_mat(:,i)]=ksdensity(Y_pred_mat(:,i));
%end

ec_pred = sum((m_pred_true>=Y_pred_CI(:,1)) & (m_pred_true<=Y_pred_CI(:,2)))/n_new;

ci_length = mean(Y_pred_CI(:,2)-Y_pred_CI(:,1));
    
%% Predictive density estimates
%[Y_fpred]=NPRegNWfprediction(Y_grid,x_new,1,1,w,theta,J);
[Y_fpred,l_fpred,u_fpred]=NPRegNWfprediction(Y_grid,x_new,q,p,w,theta,J,0.95);
%For the data
%[Y_f,l_f,u_f]=NPRegNWfprediction(Y_grid,x,q,p,w,theta,J,0.95);

% % For simulated data, we can calculate the true predictive density for
% % comparison with estimates
%Y_fpred_true=normpdf(repmat(Y_grid',n_new,1),repmat((5./(1+exp(-x_new))),1,length(Y_grid)),repmat((.25+exp((x_new-6)/3)),1,length(Y_grid)));
% For the data
%Y_f_true=normpdf(repmat(Y_grid',n,1),repmat((5./(1+exp(-x))),1,length(Y_grid)),repmat((.25+exp((x-6)/3)),1,length(Y_grid)));

% Compute L_1 error between true and estimated (mean posterior) 
% conditional density for each new covariate value 
L1_f_pred=sum(abs(Y_fpred_true'-Y_fpred),2)*(Y_grid(2)-Y_grid(1));
% Integrating with respect to the sampling distribution of the covariates
%avgL1_f_pred=normpdf(x_new',0,2.5)*L1_f_pred*(x_new(2)-x_new(1));
% Simple average
savgL1_f_pred=mean(L1_f_pred);
% Maximizing over the covariates
%maxL1_f_pred=max(L1_f_pred);
% For the data
%L1_f=sum(abs(Y_f_true-Y_f),2)*(x_new(2)-x_new(1));
% Integrating with respect to the sampling distribution of the covariates
%avgL1_f=mean(L1_f);
% Maximizing over the covariates
%maxL1_f=max(L1_f);

% %Compute 95% CI based on density estimates
% Y_pred_fCI=zeros(n_new,2);
% Y_Fpred=cumsum(Y_fpred,2)*(Y_grid(2)-Y_grid(1));
% for i=1:n_new
%     aux=1;
%     Y_Fpred(i,:)=Y_Fpred(i,:)+(1-Y_Fpred(i,ny))/2;
%     while(Y_Fpred(i,aux)<.025)
%         aux=aux+1;
%     end
%     if aux==1
%         i
%         Y_pred_fCI(i,1)=Y_grid(aux);
%     else
%         Y_pred_fCI(i,1)=(Y_grid(aux)-Y_grid(aux-1))*(0.025-Y_Fpred(i,aux-1))/(Y_Fpred(i,aux)-Y_Fpred(i,aux-1))+Y_grid(aux-1);
%     end
%     aux=ny;
%     while(Y_Fpred(i,aux)>.975)
%         aux=aux-1;
%     end
%     if aux==ny
%         i
%         Y_pred_fCI(i,2)=Y_grid(aux); 
%     else
%         Y_pred_fCI(i,2)=(Y_grid(aux+1)-Y_grid(aux))*(0.975-Y_Fpred(i,aux))/(Y_Fpred(i,aux+1)-Y_Fpred(i,aux))+Y_grid(aux);
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Model fitness via realized discrepancies %%
% % For all observations
% [Discrepancy_Y,Y_rep,Discrepancy_Y_rep,estimated_p]=NPReg_RDPredFitness(Y,x,q,p,w,theta,J);
% 
% % For "central" observations (covariate values in [-6,6] only
% x_c=x(x>-6);
% Y_c=Y(x>-6);
% Y_c=Y_c(x_c<6);
% x_c=x_c(x_c<6);
% [Discrepancy_Y_c,Y_rep_c,Discrepancy_Y_rep_c,estimated_p_c]=NPReg_RDPredFitness(Y_c,x_c,q,p,w,theta,J);
% 
% % For new data points in the complete data range
% xx=(-8:.1:10)';
% nn=size(xx,1);
% YY=5./(1+exp(-xx))+randn(nn,1).*(.25+exp((xx-6)/3));
% [Discrepancy_YY,YY_rep,Discrepancy_YY_rep,estimated_pp]=NPReg_RDPredFitness(Y_c,x_c,q,p,w,theta,J);


%% SAVE WORKSPACE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mcmc=[2000,1,2000];
save 'example1_nwreg_mmewerrors'

csvwrite('ex1_pred_nwreg_mmerrors.csv',[Y_pred,Y_pred_CI(:,1:2)]);
csvwrite('ex1_fpred_nwreg_mmerrors.csv',Y_fpred);
csvwrite('ex1_config_nwreg_mmerrors.csv',d);
csvwrite('ex1_lfpred_nwreg_mmerrors.csv',l_fpred);
csvwrite('ex1_ufpred_nwreg_mmerrors.csv',u_fpred);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 0) plot data
% figure;
% hold on;
% % Sigmoidal curve x(1)=0
% plot(x,Y, 'b.');
% 
% %with true curve
% % Sigmoidal true mean
% plot(x_new,Y_pred_true, 'b-');
% 
% %title('The data');
% xlabel('x');
% ylabel('y');
% %ylim([-2,7])
% %xlim([-6,6])
% %legend('x_1=1','x_1=0', 'Location', 'North');
% hold off
% 
% %% 1) The eight configurations with highest estimated posterior probability
% figure;
% colors=['r','g','b','m','k','c','y'];
% %markers=['o','x','*','+'];
% for j=1:1
%     subplot(1,1,1);
%     hold on;
%     uniqued=max(dR(j,:));
%     for i=1:min(uniqued,7)
%         plot(x(dR(j,:)==i),Y(dR(j,:)==i),'o','MarkerEdgeColor',colors(i));
%         %title(['Configuration post prob:',num2str(pdR(j))]);
%         xlabel('x')
%         ylabel('y')
%     end
%     if uniqued>7
%         for i=8:uniqued
%             plot(x(dR(j,:)==i),Y(dR(j,:)==i),'x','MarkerEdgeColor',colors(i-7));
%         %title(['Configuration post prob:',num2str(pdR(j))]);
%         xlabel('x')
%         ylabel('y')
%         end
%     end
% end
% % ylim([-2,7])
% % xlim([-6,6])
% 
% %% 2) Trace, autocorrelation and moving average plots of number of clusters asigned to the data 
% figure;
% uniqued=zeros(1,S);
% for i=1:S
%     uniqued(i)=size(unique(d(i,:)),2);
% end
% plot(1:S,uniqued)
% xlabel('iteration');
% ylabel('number of clusters');
% title('Trace plot of the number of clusters for the data');
% 
% nLags=100;
% figure;
% autocorr(uniqued,nLags)
% title('Autocorrelation function: number of clusters for the data');
% 
% figure;
% plot(1:S,cumsum(uniqued')./(1:S)')
% xlabel('iteration');
% ylabel('mean number of clusters');
% title('Moving average plot of the number of clusters for the data');
% 
% clear uniqued
% 
% %% 3) Trace, autocorrelation and moving average plots of the mass parameter of the Dirichlet process wen a hyperprior is used
% figure;
% plot(1:S,PP{3})
% xlabel('iteration');
% ylabel('mass parameter');
% title('Trace plot of the Dirichlet process mass parameter');
% 
% nLags=100;
% figure;
% autocorr(PP{3},nLags)
% title('Autocorrelation function: Mass parameter');
% 
% figure;
% plot(1:S,cumsum(PP{3})./(1:S)')
% xlabel('iteration');
% ylabel('Mean mass parameter');
% title('Moving average plot of the Dirichlet process mass parameter');
% 
% %% 4) Trace, autocorrelation and moving average plots of the of the precision term, tau
% figure;
% plot(1:S,theta{5})
% xlabel('iteration');
% ylabel('\tau');
% title('Trace plot of the precision term for the mixture components');
% 
% nLags=100;
% figure;
% autocorr(theta{5},nLags)
% title('Autocorrelation function: \tau');
% 
% figure;
% plot(1:S,cumsum(theta{5})./(1:S)')
% xlabel('iteration');
% ylabel('Mean \tau');
% title('Moving average plot of the precision term for the mixture components');
% 
% %% 5) Trace, autocorrelation and moving average plots of the of the regression mean for a fixed x_new value
% index=(1:S);
% index=index(((x_new(:,2)==2).*(x_new(:,1)==0))==1);
% 
% figure;
% plot(1:S,Y_pred_mat(:,index))
% xlabel('iteration');
% ylabel(['E[Y|x=(',int2str(x_new(index,1)),',',...
%     num2str(x_new(index,2)),')]']);
% title(['Trace plot of the regression mean at x=(',...
%     int2str(x_new(index,1)),',',num2str(x_new(index,2)),')']);
% 
% nLags=100;
% figure;
% autocorr(Y_pred_mat(:,index),nLags)
% title(['Autocorrelation function: E[Y|x=(',int2str(x_new(index,1)),',',...
%     num2str(x_new(index,2)),')]']);
% 
% figure;
% plot(1:S,cumsum(Y_pred_mat(:,index))./(1:S)')
% hold on
% plot([0,S],[Y_pred_true(index),Y_pred_true(index)],'r')
% hold off
% xlabel('iteration');
% ylabel(['Mean E[Y|x=(',int2str(x_new(index,1)),',',...
%     num2str(x_new(index,2)),')]']);
% title(['Moving average plot of E[Y|x=(',int2str(x_new(index,1)),',',...
%     num2str(x_new(index,2)),')]']);
% 
% %% 6) Plot the prediction with credible intervals and the true predictive mean
% % Intervals calculated directly from posterior sample 
% figure;
% hold on
% % Sigmoidal prediction x(1)=0
% plot(x_new,Y_pred,'b-');
% % Sigmoidal true mean x(1)=0
% plot(x_new,Y_pred_true, 'k-');
% % Sigmoidal lower bound x(1)=0
% plot(x_new,Y_pred_CI, 'b--');
% % Sigmoidal upper bound x(1)=0
% plot(x_new,Y_pred_CI, 'b--');
% %title('Prediction with credible intervals and true mean')
% xlabel('x')
% ylabel('y')
% hold off
% 
% %% 7) Plot the prediction with credible intervals and the true predictive mean
% % Intervals calculated from the estimated predictive densities
% figure;
% plot(x_new(x_new(:,1)==0,2),Y_pred(x_new(:,1)==0),'r-');
% hold on
% plot(x_new(x_new(:,1)==1,2),Y_pred(x_new(:,1)==1),'b-');
% %plot(x(:,2),Y, 'rx');
% plot(x_new(x_new(:,1)==0,2),Y_pred_true(x_new(:,1)==0),'k-');
% plot(x_new(x_new(:,1)==1,2),Y_pred_true(x_new(:,1)==1),'k-');
% plot(x_new(x_new(:,1)==0,2),Y_pred_fCI(x_new(:,1)==0,1),'r--');
% plot(x_new(x_new(:,1)==0,2),Y_pred_fCI(x_new(:,1)==0,2),'r--');
% plot(x_new(x_new(:,1)==1,2),Y_pred_fCI(x_new(:,1)==1,1),'b--');
% plot(x_new(x_new(:,1)==1,2),Y_pred_fCI(x_new(:,1)==1,2),'b--');
% %title('Data and prediction')
% %legend('Pred. (x_1=0)','Pred. (x_1=1)','Data','True mean', 'Location','North')
% xlabel('x_2')
% ylabel('y')
% hold off
% 
% %% 8) Estimated conditional density surface
% load('MyColormaps','inverse_hot_hsv')
% 
% x1toplot=1;
% figure;
% colormap(inverse_hot_hsv)
% surf(x_new(x_new(:,1)==x1toplot,2)-0.05,Y_grid-0.05,...
%     (Y_fpred(x_new(:,1)==x1toplot,:))')
% axis([-5 5 -10 10 -1 0.5])
% hold on
% surf(x_new(x_new(:,1)==x1toplot,2)-0.05,Y_grid-0.05,...
%     -1*ones(size((Y_fpred(x_new(:,1)==x1toplot,:))')),...
%     (Y_fpred(x_new(:,1)==x1toplot,:))')
% shading interp
% hold off
% xlabel('x')
% ylabel('y')
% zlabel('f(y|x)')
% title(['Estimated conditional density surface for x_1=',int2str(x1toplot)])
% colorbar('location','eastoutside')
% %caxis(caxislim)
% 
% %% 9) Color plot of the estimated posterior density of the regression curve with mean and true
% load('MyColormaps','inverse_hot_hsv')
% 
% figure;
% pcolor(x_new,Y_mat,Y_pred_hist); 
% shading flat
% colormap(inverse_hot_hsv)
% %caxislim=caxis;
% hold on
% % Estimated
% plot(x_new,Y_pred,'b','LineWidth',2);
% % True
% plot(x_new,Y_pred_true,'k','LineWidth',2);
% hold off
% title('Posterior density of the regression mean');
% xlabel('x')
% ylabel('y')
% colorbar('location','eastoutside')
% % inverse_hot_hsv=colormap;
% % save('MyColormaps','inverse_hot_hsv')
% ylim([-2,8])
% 
% %% 10) Color plot of the estimated predictive density with estimated and true regression curve
% load('MyColormaps','inverse_hot_hsv')
% 
% figure;
% fig=pcolor(x_new-0.05,Y_grid-0.05,Y_fpred'); 
% shading flat
% colormap(inverse_hot_hsv)
% %colormap(hot)
% %caxislim=caxis;
% hold on
% % Estimated
% plot(x_new,Y_pred,'b','LineWidth',2);
% % True
% %plot(x_new,Y_pred_true,'k','LineWidth',2);
% %Data
% plot(x,Y, 'b.');
% hold off
% %title('Estimated conditional densities');
% xlabel('x')
% ylabel('y')
% colorbar('location','eastoutside')
% % inverse_hot_hsv=colormap;
% % save('MyColormaps','inverse_hot_hsv')
% 
% %% 10b) Color plot of the true predictive density with estimated and true regression curve
% load('MyColormaps','inverse_hot_hsv')
% 
% figure;
% pcolor(x_new-0.05,Y_grid-0.05,Y_fpred_true'); 
% caxis=caxislim;
% shading flat
% colormap(inverse_hot_hsv)
% %colormap(hot)
% %caxislim=caxis;
% hold on
% %Data
% plot(x,Y, 'b.');
% % True
% plot(x_new,Y_pred_true,'k','LineWidth',2);
% hold off
% %title('Estimated conditional densities');
% xlabel('x')
% ylabel('y')
% colorbar('location','eastoutside')
% % inverse_hot_hsv=colormap;
% % save('MyColormaps','inverse_hot_hsv')
% 
% %% 11) Plot true predictive density and estimates for some values of x
% figure;
% x_toplot=[21:20:101];
% hold on
% for j=1:5
%     %subplot(2,5,j);
%     %hold on;
%     plot(Y_grid,Y_fpred_true(x_toplot(j),:)','k');
%     plot(Y_grid,Y_fpred(x_toplot(j),:)','b');
%     title(['x_1=',int2str(x_new(x_toplot(j),1)),'; x_2=',num2str(x_new(x_toplot(j),2))])
%     xlabel('y')
%     ylabel('p(y|x)')
% end
% 
% % with credible intervals
% figure;
% x_toplot=[1:25:101,102:25:202];
% for j=1:10
%     subplot(2,5,j);
%     hold on;
%     plot(Y_grid,Y_fpred_true(x_toplot(j),:)','k');
%     plot(Y_grid,Y_fpred(x_toplot(j),:)','b');
%     plot(Y_grid,l_fpred(x_toplot(j),:)','b--');
%     plot(Y_grid,u_fpred(x_toplot(j),:)','b--');
%     title(['x_1=',int2str(x_new(x_toplot(j),1)),'; x_2=',num2str(x_new(x_toplot(j),2))])
%     xlabel('y')
%     ylabel('p(y|x)')
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% 12) Examing the covariate dependent weights for the last iteration
% it=last_s_dR(1);
% 
% % x_grid=[(-5:.1:5)';(-5:.1:5)'];
% % n_grid=size(x_grid,1);
% % x_grid=[[zeros(n_grid/2,1);ones(n_grid/2,1)],x_grid];
% 
% x_grid=x_new;
% n_grid=n_new;
%  %Compute matrix of weights for each x
% lg2=zeros(n_grid,J(it));
% if p>0
%     for h=1:p
%         lg2=lg2-0.5*theta{4}(it,h)*(kron(ones(1,J(it)),x_grid(:,q+h))-kron(ones(n_grid,1),theta{3}{it}(h,:))).^2;
%     end
% end
% wx=kron(ones(n_grid,1),w{it}).*exp(lg2); %weights *normal kernel for x_new in component j
% wx=wx./kron(ones(1,J(it)),sum(wx,2)); %normalize
% 
% % Plot first 7 weights as a function of x(2) for x(1)=1
% figure;
% hold on;
% for weight_num=1:7
%     plot(x_new,wx(:,weight_num), colors(weight_num));
% end
% 
% xlabel('x');
% ylabel('w_j(x)');
% 
% %% 13) Plot density estimates for some values of x
% colors=['m','r','g','b','k'];
% x_toplot=[21:20:101];
% figure;
% hold on;
% for i=1:5
%     plot(Y_grid,Y_fpred(x_toplot(i),:)','Color',colors(i));
%     plot(Y_grid,Y_fpred_true(x_toplot(i),:)','--','Color',colors(i));
% end
% %axis([3 12 0 0.6])
% xlabel('y')
% ylabel('f(y|x)')
% legend('x=-4','x=-2','x=0','x=2','x=4')
% hold off;
% 
% %% 14) Scatterplot of predicted vs. realized discrepancies
% % For complete data
% figure; hold on
% plot(Discrepancy_YY,Discrepancy_YY_rep,'.')
% plot([0,200],[0,200],'r')
% 
% % For "central" observations (covariate values in [-6,6] only
% figure; hold on
% plot(Discrepancy_Y_c,Discrepancy_Y_rep_c,'.')
% plot([0,200],[0,200],'r')