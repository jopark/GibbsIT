rng('shuffle')
p=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99]; % confidence levels

%% for E. coli and Mammalian cells
% 5 GAPD:PGK
% 7 PGM:ENO
% 10 PGI:PFK
% 11 PGI:ENO

%% for Clostridia
% 5 GAPD:PGK
% 7 PGM:ENO
% 11 PGI:ENO

%% Choose from above
rxn1=10;

%% for E. coli and Mammalian cells
% 12 FBA
% 13 TPI

%% for Clostridia
% 7 PGI
% 8 PFK
% 9 FBA
% 10 TPI

%% Choose from above
rxn2=13;

%% Script for boostrapping starts
% dg1 corresponds to dG of N-, P- limitation, anaerobic E. coli, untreated
% iBMK mammalian cells, and Clostridium acetobutylicum.
% dg2 corresponds to dG of N-, P-, O- upshift E. coli, oligomycin treated
% iBMK mammalian cells, and Clostridium cellulolyticum.

mindg1=-30; % the lowerbound for the pre-perturbation rxn dG (should be lowerbound of 99% confidence interval)
mindg2=-30; % the lowerbound for the post-perturbation rxn dG (should be lowerbound of 99% confidence interval)
mindg=min(mindg1,mindg2);
step=-mindg/100;

samples=1000;
dg1=zeros(1,samples);
dg2=zeros(1,samples);

for i=1:samples

%   Use the following two for lumped reactions with index: rxn1
    dg1_=sampledg(lub1,rxn1,-200*mindg,p,mindg1,0);
    dg2_=sampledg(lub2,rxn1,-200*mindg,p,mindg2,0);

%   Use the following two for individual reactions with index: rxn2
%     dg1_=sampledg2(lb1,ub1,rxn2,-200*mindg,p,mindg1,0);
%     dg2_=sampledg2(lb2,ub2,rxn2,-200*mindg,p,mindg2,0);

    dg1(i)=mean(dg1_(1:3));
    dg2(i)=mean(dg2_(1:3));
end

% Plot histogram showing the two distributions
map = brewermap(2,'Set1');
figure
[n1,x1]=histf(dg1,mindg-step/2:step:0,'facecolor',map(1,:),'facealpha',.5,'edgecolor',map(1,:));
hold on
[n2,x2]=histf(dg2,mindg-step/2:step:0,'facecolor',map(2,:),'facealpha',.5,'edgecolor',map(2,:));
box off
hold off

%% Find area overlap
overlap=sum(min(n1,n2))/samples;