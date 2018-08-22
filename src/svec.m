function [obs,obs_i,obs_l,meas,cur,svar] = svec(mea,sim,var)
%% create a single vector of all observed labeling patterns
% Input: mea (structure of measurements)
% 
% Output: obs (vector of all observed labeling fractions),obs_i and obs_l
% (indices vector and length required to decipher obs)
meas=fieldnames(mea);
meas_l=length(meas);
obs_i=zeros(meas_l+1,1);
sobs_i=zeros(meas_l+1,1);
for i=1:meas_l
	if nargin<2 || isempty(sim)
		len=length(mea.(meas{i}));
	else
		len=length(sim.(meas{i}));
		slen=length(mea.(meas{i}));
		sobs_i(i+1)=sobs_i(i)+slen;
	end
    obs_i(i+1)=obs_i(i)+len;
end
obs_l=obs_i(end);
obs=zeros(obs_l,1);
cur=zeros(obs_l,1);
if isscalar(var)
	svar=diag(var*ones(max(obs_l,sobs_i(end)),1));
else
	svar=var;
end
if nargin>1 && ~isempty(sim) && ~isempty(var)
	for i=meas_l:-1:1
		lenm=length(mea.(meas{i}));
		lens=length(sim.(meas{i}));
		svar(sobs_i(i)+1+lens:sobs_i(i)+lenm,:)=[];
		svar(:,sobs_i(i)+1+lens:sobs_i(i)+lenm)=[];
	end
end

for i=1:meas_l
	if nargin<2 || isempty(sim)
		obs(obs_i(i)+1:obs_i(i+1))=mea.(meas{i})(:);
		cur=[];
	else
		obs(obs_i(i)+1:obs_i(i+1))=mea.(meas{i})(1:obs_i(i+1)-obs_i(i))';
		cur(obs_i(i)+1:obs_i(i+1))=sim.(meas{i})(:);
	end
end