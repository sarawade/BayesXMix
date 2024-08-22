function [CS,dR,pdR,D_indexS]=relabel(d)
%  function [CS,dR,pdR,D_indexS]=relabel(d)
%
% Relabels vector of component indices so that the Y(1) is
% always allocated to the first component, i.e. dR(s,1)=1 for all s.
% The first Y(i) with i>1 such that d(s,i) is not equal to d(s,1) is allocated
% to the second component, i.e. dR(i,s)=2, and so on...
% Then configurations are sorted in descending order according to their
% estimated posterior probabilities and redundant configurations are
% eliminated
% INPUT:  
%   D: S x n    matrix of indices for each data point sampled from
%       the posterior
%   
% OUTPUT: 
%   CS: Number of unique Configurations among the S posterior samples
%       (CS<=S)
%   dR: CS x n matrix unique configurations. Each 1 x n vector is a
%       configuration relabeled by order of appearance. The vectors are then
%       sorted by their estimated posterior probabilities.
%   pdR: CS x 1 vector of estimated posterior probability of each configuration in dR
%       pdR(i)=frequency of dR(i,:)/S
%   D_indexS: indices at which the unique configurations were last visited,
%       ordered by configuration frequency
  

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE AND ALLOCATE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[S,n]=size(d);
D_count=ones(S,1);
D_index=S;
CS=1;
D_ordered=zeros(S,n);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REORDER LABELS AND CALCULATE CONFIGURATION FREQUENCY %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for s=S:-1:1
	% Reorder the configuration
	uniq_s=unique(d(s,:));
	ks=length(uniq_s);
	for h=1:ks
		D_ordered(s,d(s,:)==uniq_s(h))=h;
	end
	if s<S
		% Check if we have seen this configuration before (for a smaller s)
		before=sum(D_ordered(D_index,:)==kron(ones(length(D_index),1),D_ordered(s,:)),2)==n;
		if sum(before)>0 % We have seen it before
			D_count(before)=D_count(before)+1;
        else % We have not seen it before
			CS=CS+1;
			D_index(CS)=s;
		end
	end	
end
D_count=D_count(1:CS);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SORT CONFIGURATIONS BY FREQUENCY AND ELIMINATE REDUNDANCY %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[D_countS,ID_countS]=sort(D_count,'descend');
D_indexS=D_index(ID_countS);
dR=D_ordered(D_indexS,:);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE POSTERIOR PROBABILITY ESTIMATES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pdR=D_countS/S;

%it=D_indexS(1);
