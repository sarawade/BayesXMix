function mass=massupdate(massa,b,n,da,ka,Da)
% function mass=massupdate(massa,b,n,da,ka,Da)
% Updates the mass parameter of the Dirichlet process prior with Gamma
% hyperprior
% INPUT:
% massa: current value for the mass parameter
% b: 2x1 vector of hyperparameters for the Gamma hyperprior
% n: sample size
% da: 1 x n vector with current state of index variables
% ka: 1 x n vector with current k values
% Da: cell array of length n with current index values for latent variables
%
% OUTPUT:
% wa: vector of length Ja of updated mixture weights
% mass: updated mass parameter for the Dirichlet process

%% Find the number of unique components currently assigned to data points or latent variables
Ds=unique(da);
for i=1:n
    if ka(i)>0
        Ds=unique([Ds,Da{i}]);
    end
end
nDs=size(Ds,2);
N=n+sum(ka);

%% Sample auxiliary variable
e=betarnd(massa+1,N);

%% Updated hyperparameters
z1=b(1)+nDs-1;
z2=b(2)-log(e);

%% Sample new mass parameter from a mixture of Gamma densities
r=z1/(z1+N*z2);
if rand<r
    mass=gamrnd(z1+1,1/z2);
else
    mass=gamrnd(z1,1/z2);
end

if mass<=0.1
    mass
end