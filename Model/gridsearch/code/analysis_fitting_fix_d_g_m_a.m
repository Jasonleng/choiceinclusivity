clear all; close all;

files = dir('../outcome/*.mat');


%%
tic;
% datalist = {};
% 
% for i = 1:length(files)
%     datalist{i} = load(['../outcome/' files(i).name]);
%     
%     data.I_RT(i,:,:,:,:) = datalist{i}.I_RT;
%     data.I_Acc(i,:,:,:,:) = datalist{i}.I_Acc;
%     data.b1_RT(i,:,:,:,:) = datalist{i}.b1_RT;
%     data.b2_RT(i,:,:,:,:) = datalist{i}.b2_RT;
%     data.b1_Acc(i,:,:,:,:) = datalist{i}.b1_Acc;
%     data.b2_Acc(i,:,:,:,:) = datalist{i}.b2_Acc;
%     data.decay(i,:,:,:,:) = datalist{i}.decay;
%     data.mutualInhib(i,:,:,:,:) = datalist{i}.mutualInhib;
%     data.thresh(i,:,:,:,:) = datalist{i}.thresh;
%     data.gain(i,:,:,:,:) = datalist{i}.gain;
%     data.nDT(i) = datalist{i}.nDT;
%     data.order(i) = datalist{i}.order;
% end
% datalist = {};

for i = 1:length(files)
    datalist = load(['../outcome/' files(i).name]);
    
    data.I_RT(i,:,:,:,:) = datalist.I_RT;
    data.I_Acc(i,:,:,:,:) = datalist.I_Acc;
    data.b1_RT(i,:,:,:,:) = datalist.b1_RT;
    data.b2_RT(i,:,:,:,:) = datalist.b2_RT;
    data.b1_Acc(i,:,:,:,:) = datalist.b1_Acc;
    data.b2_Acc(i,:,:,:,:) = datalist.b2_Acc;
    data.decay(i,:,:,:,:) = datalist.decay;
    data.mutualInhib(i,:,:,:,:) = datalist.mutualInhib;
    data.thresh(i,:,:,:,:) = datalist.thresh;
    data.gain(i,:,:,:,:) = datalist.gain;
    data.nDT(i) = datalist.nDT;
    data.order(i) = datalist.order;
end


toc;
%%


addpath('../../models/');

coeff_rt = readtable('../../data/outcome/coeff_rt_12.csv');
coeff_acc = readtable('../../data/outcome/coeff_acc_12.csv');
intercept = readtable('../../data/outcome/intercept_acc_rt_12.csv');

study = 'Study 12';

% 
% c = [];
% order = [];
% nDT = [];
% file = {};

% % maxCoeffRT = [];
% maxCoeffAcc = [];
% maxRT = [];
% maxAcc = [];

decays = [0 0];
mis = [0 0];
threshs = [0 0];
gains = [0 0];
orders = [0 0];
nDTs = [0 0];
nondecsT = [0 0];
coeff_rt_sub_or = coeff_rt(strcmp(coeff_rt.study,study) & strcmp(coeff_rt.catChoiceCondition,'OR'),:);
coeff_acc_sub_or = coeff_acc(strcmp(coeff_acc.study,study) & strcmp(coeff_acc.catChoiceCondition,'OR'),:);
intercept_sub_or = intercept(strcmp(intercept.study,study) & strcmp(intercept.catChoiceCondition,'OR'),:);
coeff_rt_sub_xor = coeff_rt(strcmp(coeff_rt.study,study) & strcmp(coeff_rt.catChoiceCondition,'XOR'),:);
coeff_acc_sub_xor = coeff_acc(strcmp(coeff_acc.study,study) & strcmp(coeff_acc.catChoiceCondition,'XOR'),:);
intercept_sub_xor = intercept(strcmp(intercept.study,study) & strcmp(intercept.catChoiceCondition,'XOR'),:);

%% Fit XOR

nDTlist = 60:10:160;

mincost = 100000;
for order = 1:0.1:2
    
    for m = 1:11
        for k = 1:11

            curDataMask1 = data.order == order & data.nDT == nDTlist(m);
            curDataMask2 = data.order == order & data.nDT == nDTlist(k);

            
            costmat = zeros(11,11,11,11,12);

            bestndt = -((squeeze(data.I_RT(curDataMask1,:,:,:,:)) - intercept_sub_xor.I_rt)/intercept_sub_xor.se_rt^2 + ...
                (squeeze(data.I_RT(curDataMask2,:,:,:,:)) - intercept_sub_or.I_rt)/intercept_sub_or.se_rt^2)/(1/intercept_sub_xor.se_rt^2 + ...
                1/intercept_sub_or.se_rt^2);
            
            bestndt(bestndt<0) = 0;
            costmat(:,:,:,:,1) = ((squeeze(data.b1_RT(curDataMask1,:,:,:,:)) - coeff_rt_sub_xor.beta_OV)/coeff_rt_sub_xor.se_OV).^2;
            costmat(:,:,:,:,2) = ((squeeze(data.b2_RT(curDataMask1,:,:,:,:)) - coeff_rt_sub_xor.beta_VD)/coeff_rt_sub_xor.se_VD).^2;
            costmat(:,:,:,:,3) = ((squeeze(data.b1_Acc(curDataMask1,:,:,:,:)) - coeff_acc_sub_xor.beta_OV)/coeff_acc_sub_xor.se_OV).^2;
            costmat(:,:,:,:,4) = ((squeeze(data.b2_Acc(curDataMask1,:,:,:,:)) - coeff_acc_sub_xor.beta_VD)/coeff_acc_sub_xor.se_VD).^2;
            costmat(:,:,:,:,5) = ((bestndt + squeeze(data.I_RT(curDataMask1,:,:,:,:)) - intercept_sub_xor.I_rt)/intercept_sub_xor.se_rt).^2;
            costmat(:,:,:,:,6) = ((squeeze(data.I_Acc(curDataMask1,:,:,:,:)) - intercept_sub_xor.I_acc)/intercept_sub_xor.se_acc).^2;
            costmat(:,:,:,:,7) = ((squeeze(data.b1_RT(curDataMask2,:,:,:,:)) - coeff_rt_sub_or.beta_OV)/coeff_rt_sub_or.se_OV).^2;
            costmat(:,:,:,:,8) = ((squeeze(data.b2_RT(curDataMask2,:,:,:,:)) - coeff_rt_sub_or.beta_VD)/coeff_rt_sub_or.se_VD).^2;
            costmat(:,:,:,:,9) = ((squeeze(data.b1_Acc(curDataMask2,:,:,:,:)) - coeff_acc_sub_or.beta_OV)/coeff_acc_sub_or.se_OV).^2;
            costmat(:,:,:,:,10) = ((squeeze(data.b2_Acc(curDataMask2,:,:,:,:)) - coeff_acc_sub_or.beta_VD)/coeff_acc_sub_or.se_VD).^2;
            costmat(:,:,:,:,11) = ((bestndt + squeeze(data.I_RT(curDataMask2,:,:,:,:)) - intercept_sub_or.I_rt)/intercept_sub_or.se_rt).^2;
            costmat(:,:,:,:,12) = ((squeeze(data.I_Acc(curDataMask2,:,:,:,:)) - intercept_sub_or.I_acc)/intercept_sub_or.se_acc).^2;
            bestndtmat(:,:,:,:) = bestndt;
            
            cost = sum(costmat,5);
            % cost = max(costmat,[],5);
            if min(cost(:)) < mincost
                [v,loc] = min(cost(:));
                mincost = v;
                [ii,jj,kk,ww] = ind2sub(size(cost),loc);
                for j = 1:2
                    decays(j) = data.decay(1,ii,jj,kk,ww);
                    threshs(j) = data.thresh(1,ii,jj,kk,ww);
                    gains(j) = data.gain(1,ii,jj,kk,ww);
                    mis(j) = data.mutualInhib(1,ii,jj,kk,ww);
                    orders(j) = order;
                    nondecsT(j) = bestndtmat(ii,jj,kk,ww);
                end
                nDTs(1) = nDTlist(m);
                nDTs(2) = nDTlist(k);
            end
        end
    end
end
%%

% besttheta = [bestgains(1) 0 bestdecays(1) bestdecays(1)*bestmis(1) bestthreshs(1) bestnondecsT(1) 0.01 0.5;
%             bestgains(2) 0 bestdecays(2) bestdecays(2)*bestmis(2) bestthreshs(2) bestnondecsT(2) 0.01 0.5];


besttheta = [gains(1) 0 decays(1) decays(1)*mis(1) threshs(1) nondecsT(1) 0.01 0.5;
             gains(2) 0 decays(2) decays(2)*mis(2) threshs(2) nondecsT(2) 0.01 0.5];



outcome = summaryPredictionSeparate(besttheta,nDTs,orders);

%%
% file = file';
% 
% summary = table(file,c,order,nDT,maxCoeffAcc,maxRT,maxAcc);
% 
% summary.product = log(summary.maxAcc) + log(summary.maxRT);

% %%
% 
% data = load('fitting_LCA_EXT/5396632/LCA_RT_quad_5396632_task_20.mat'); % Best selected one
% 
% aggregation = cat(3,data.b1_Acc_mat,data.b1_RT_mat,data.b2_Acc_mat,data.b2_RT_mat);
% 
% sumMat = max(aggregation,[],3);
% 
% sumMat == min(sumMat(:))
% 
% %%
% params = [data.decay_mat(sumMat == min(sumMat(:))),...
%     0,data.thresh_mat(sumMat == min(sumMat(:))),data.c,0.2,0.05,0.5];
% outcome = summaryPrediction([data.decay_mat(sumMat == min(sumMat(:))),...
%     0,data.thresh_mat(sumMat == min(sumMat(:))),data.c,0.2,0.05,0.5],...
%     data.nDT,data.order,[0,1]);
% 
% outcome
% %%
table_outcome = array2table(outcome');

table_outcome.Properties.VariableNames(1:6) = {'I_RT','b1_RT','b2_RT','I_Acc','b1_Acc','b2_Acc'};

table_outcome.condition = {'XOR';'OR'};

writetable(table_outcome,'../result/prediction_task_12_ddl.csv');