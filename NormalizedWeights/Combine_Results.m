%% For combining results from consecutive runs

%% Load results from first file
clear
WS=load('example3_n200_S10000_b0_ini-test_hp_LS');

% Extract data
Y=WS.Y;
x=WS.x;

% Extract parameters
mcmc=WS.mcmc;
phi=WS.phi;
PP=WS.PP;
initialState=WS.initialState;
lastState=WS.lastState;
cputime=WS.cputime;
Hyperparameters=WS.Hyperparameters;

% Extract results
theta=WS.theta;
d=WS.d;
w=WS.w;
J=WS.J;

% Extract preprocessing results
x_new=WS.x_new;
n_new=WS.n_new;
Y_pred_mat=WS.Y_pred_mat;
Y_pred_true=WS.Y_pred_true;
Y_grid=WS.Y_grid;
ny=WS.ny;
Y_fpred_true=WS.Y_fpred_true;

% Delete what is not needed
clear WS

%% Load results from second file and add to previous
WS=load('example3_n200_S10000_b10000_ini-test_hp_LS');
mcmc(3)=mcmc(3)+WS.mcmc(3);
d=[d;WS.d];
w=[w;WS.w];
J=[J;WS.J];
for i=1:size(theta,1)
    theta{i,1}=[theta{i,1};WS.theta{i,1}];
end
if strcmp('DirichletHP',PP{1})
    PP{3}=[PP{3},WS.PP{3}];
end

% Load preprocessing results and add to previous
Y_pred_mat=[Y_pred_mat;WS.Y_pred_mat;];

% Delete what is not needed
cputime=cputime+WS.cputime;
clear WS i

%% Cutting to keep only last iterations (for longer burn-in)
S=1000; 
mcmc(1)=mcmc(1)+(mcmc(3)-S)*mcmc(2);
mcmc(3)=S;
for i=1:size(theta,1)
    theta{i,1}=theta{i,1}(end-S+1:end);
end
d=d(end-S+1:end,:);
w=w(end-S+1:end);
J=J(end-S+1:end);
Y_pred_mat=Y_pred_mat(end-S+1:end,:);
if strcmp('DirichletHP',PP{1})
    PP{3}=PP{3}(end-S+1:end);
end
clear i 

%% Thinning
mcmc(2)=20; %thinning
iterations=1:mcmc(2):mcmc(3);
S=size(iterations,2);
for i=1:size(theta,1)
    theta{i,1}=theta{i,1}(iterations);
end
d=d(iterations,:);
w=w(iterations);
J=J(iterations);
Y_pred_mat=Y_pred_mat(iterations,:);
if strcmp('DirichletHP',PP{1})
    PP{3}=PP{3}(iterations');
end
clear i iterations

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recalculate (post process)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clustering
[CS,dR,pdR,last_s_dR]=relabel(d);

%% Regression curve
% mean
Y_pred=(sum(Y_pred_mat,1)/S)';

% confidence intervals
cred_prob=0.95;
predS_sort=sort(Y_pred_mat); 
l_pred=predS_sort(uint16(S*(1-cred_prob)/2),:); % Lower bound
u_pred=predS_sort(uint16(S*(cred_prob+(1-cred_prob)/2)),:); % Upper bound
Y_pred_CI=[l_pred',u_pred'];
clear predS_sort u_pred l_pred

% estimated posterior density (smoothed histograms)
Y_pred_hist=zeros(100,n_new);
Y_mat=zeros(100,n_new);
for i=1:n_new
    [Y_pred_hist(:,i),Y_mat(:,i)]=ksdensity(Y_pred_mat(:,i));
end

%% Predictive densities

%Densities with credible bands
[Y_fpred,l_fpred,u_fpred]=NPRegNWfprediction(Y_grid,x_new,1,1,w,theta,J,0.95);

%Compute 95% CI based on density estimates
Y_pred_fCI=zeros(n_new,2);
Y_Fpred=cumsum(Y_fpred,2)*(Y_grid(2)-Y_grid(1));
for i=1:n_new
    aux=1;
    Y_Fpred(i,:)=Y_Fpred(i,:)+(1-Y_Fpred(i,ny))/2;
    while(Y_Fpred(i,aux)<.025)
        aux=aux+1;
    end
    if aux==1
        i
        Y_pred_fCI(i,1)=Y_grid(aux);
    else
        Y_pred_fCI(i,1)=(Y_grid(aux)-Y_grid(aux-1))*(0.025-Y_Fpred(i,aux-1))/(Y_Fpred(i,aux)-Y_Fpred(i,aux-1))+Y_grid(aux-1);
    end
    aux=ny;
    while(Y_Fpred(i,aux)>.975)
        aux=aux-1;
    end
    if aux==ny
        i
        Y_pred_fCI(i,2)=Y_grid(aux); 
    else
        Y_pred_fCI(i,2)=(Y_grid(aux+1)-Y_grid(aux))*(0.975-Y_Fpred(i,aux))/(Y_Fpred(i,aux+1)-Y_Fpred(i,aux))+Y_grid(aux);
    end
end

%% Save workspace
save 'example3_n200_S15000_b25000_th20_hp'