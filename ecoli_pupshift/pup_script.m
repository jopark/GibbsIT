%% flux and confidence interval estimation script

load ec_pup
[mea,covar,raw,avgvar]=xlsreadplus(xlsname,{'5D','1_2_13C'});

% account for low variance to avoid overfitting
var=diag(covar);
p25=prctile(var(var~=0),25);
p50=prctile(var(var~=0),50);
p75=prctile(var(var~=0),75);
minvar=0.004^2;
maxvar=0.035^2';
var(var==0)=0.01^2;
var(var<p25)=p25;
var(var<minvar)=minvar;
var(var>p75)=p75;
var(var>maxvar)=maxvar;
var=diag(var);

%%
% input metabolites
input.GLYCOGEN_IN1__12C=isotopomervector([1 1 0 0 0 0 0 0 0 0 0 0 0],1);
input.GLYCOGEN_IN2__12C=isotopomervector([0 0 0 0 0 0 0 0 0 0 0 0 0],1);
input.GLYCOGEN_IN1__5H=isotopomervector([0 0 0 0 0 0 0 0 0 0 0 0 0],1);
input.GLYCOGEN_IN2__5H=isotopomervector([0 0 0 0 0 0 0 0 0 0 0 0 0],1);
input.GLC_IN__12C=isotopomervector([1 1 0 0 0 0 0 0 0 0 0 0 0],1);
input.GLC_IN__5H=isotopomervector([0 0 0 0 0 0 0 0 0 0 1 0 0],1);
input.CO2_IN__12C=isotopomervector([0],1);
input.CO2_IN__5H=isotopomervector([0],1);
input.NADH_IN__12C=isotopomervector([0],1);
input.NADH_IN__5H=isotopomervector([0],1);
input.NADPH_IN__12C=isotopomervector([0],1);
input.NADPH_IN__5H=isotopomervector([0],1);
input.H_IN__12C=isotopomervector([0],1);
input.H_IN__5H=isotopomervector([0],1);
input.PYR_IN__12C=isotopomervector([0 0 0 0 0 0],1);
input.PYR_IN__5H=isotopomervector([0 0 0 0 0 0],1);
input.RIBOSE_IN__12C=isotopomervector([0 0 0 0 0 0 0 0 0 0 0],1);
input.RIBOSE_IN__5H=isotopomervector([0 0 0 0 0 0 0 0 0 0 0],1);

free_net=rand(size(model.kernel_net,2),1);
free_xch=rand(size(model.kernel_xch,2),1);

[fmea,fmeaStr]=xlsreadfmea(xlsname,model,free_net,free_xch);
[ineq,ineqStr]=xlsreadineq(xlsname,model,free_net,free_xch);

simulate=str2func(xmlfile);

%% evaluate the best fit fluxes and labeling pattern
[net_opt,xch_opt,info,net_,xch_,info_,fval_,I]=pargloptflux(simulate,model,free_net,free_xch,ineq,[],input,mea,fmea,var,100);

[tflux,tiso,tscore]=writetableplus(simulate,model,net_opt,xch_opt,input,mea,fmea,var);
sim=simulate(net_opt,xch_opt,input);

%% estimate confidence interval
[lb,ub,hs,net_nopt,xch_nopt,nscore,nflag,I,lub_lg,lub_lg_s]=parconfest3_1lg2(simulate,model,net_opt,xch_opt,ineq,input,mea,fmea,var);

lub_lgg=[lub_lg_s(:) num2cell([-8.314*310.15/1000*log(lub_lg(:,2)) -8.314*310.15/1000*log(lub_lg(:,1)) lub_lg])];
ub=[ub -8.314*310.15/1000*log(lb(:,5))];
lb=[lb -8.314*310.15/1000*log(ub(:,5))];
