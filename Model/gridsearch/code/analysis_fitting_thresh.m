clear all; close all;

files = dir('../outcome/*.mat');

coeff_rt = readtable('../../data/outcome/coeff_rt_12.csv');
coeff_acc = readtable('../../data/outcome/coeff_acc_12.csv');
intercept = readtable('../../data/outcome/intercept_acc_rt_12.csv');

study = 'Study 12';

% 
% c = [];
% order = [];
% nDT = [];
% file = {};

% maxCoeffRT = [];
maxCoeffAcc = [];
maxRT = [];
maxAcc = [];
mincost = 100000;
decays = [0 0];
mis = [0 0];
threshs = [0 0];
gains = [0 0];
orders = [0 0];
nDTs = [0 0];
nondecsT = [0 0];
file = '';
coeff_rt_sub_or = coeff_rt(strcmp(coeff_rt.study,study) & strcmp(coeff_rt.catChoiceCondition,'OR'),:);
coeff_acc_sub_or = coeff_acc(strcmp(coeff_acc.study,study) & strcmp(coeff_acc.catChoiceCondition,'OR'),:);
intercept_sub_or = intercept(strcmp(intercept.study,study) & strcmp(intercept.catChoiceCondition,'OR'),:);
coeff_rt_sub_xor = coeff_rt(strcmp(coeff_rt.study,study) & strcmp(coeff_rt.catChoiceCondition,'XOR'),:);
coeff_acc_sub_xor = coeff_acc(strcmp(coeff_acc.study,study) & strcmp(coeff_acc.catChoiceCondition,'XOR'),:);
intercept_sub_xor = intercept(strcmp(intercept.study,study) & strcmp(intercept.catChoiceCondition,'XOR'),:);
for i = 1:length(files)
    data = load(['../outcome/' files(i).name]);

    costmat = zeros(11,11,11,11,11,12);
    bestndtmat = zeros(11,11,11,11,11);
    for j = 1:11
        for k = 1:11
%         
%         rs = [1 3];
%         
%             coeff_rt_sub = coeff_rt(strcmp(coeff_rt.study,study) & strcmp(coeff_rt.catChoiceCondition,conds{j}),:);
%             coeff_acc_sub = coeff_acc(strcmp(coeff_acc.study,study) & strcmp(coeff_acc.catChoiceCondition,conds{j}),:);
%             intercept_sub = intercept(strcmp(intercept.study,study) & strcmp(intercept.catChoiceCondition,conds{j}),:);
            
            bestndt = -((squeeze(data.I_RT(:,:,j,:)) - intercept_sub_xor.I_rt)/intercept_sub_xor.se_rt^2 + ...
                (squeeze(data.I_RT(:,:,k,:)) - intercept_sub_or.I_rt)/intercept_sub_or.se_rt^2)/(1/intercept_sub_xor.se_rt^2 + ...
                1/intercept_sub_or.se_rt^2);
            costmat(:,:,:,j,k,1) = ((squeeze(data.b1_RT(:,:,j,:)) - coeff_rt_sub_xor.beta_OV)/coeff_rt_sub_xor.se_OV).^2;
            costmat(:,:,:,j,k,2) = ((squeeze(data.b2_RT(:,:,j,:)) - coeff_rt_sub_xor.beta_VD)/coeff_rt_sub_xor.se_VD).^2;
            costmat(:,:,:,j,k,3) = ((squeeze(data.b1_Acc(:,:,j,:)) - coeff_acc_sub_xor.beta_OV)/coeff_acc_sub_xor.se_OV).^2;
            costmat(:,:,:,j,k,4) = ((squeeze(data.b2_Acc(:,:,j,:)) - coeff_acc_sub_xor.beta_VD)/coeff_acc_sub_xor.se_VD).^2;
            costmat(:,:,:,j,k,5) = ((bestndt + squeeze(data.I_RT(:,:,j,:)) - intercept_sub_xor.I_rt)/intercept_sub_xor.se_rt).^2;
            costmat(:,:,:,j,k,6) = ((squeeze(data.I_Acc(:,:,j,:)) - intercept_sub_xor.I_acc)/intercept_sub_xor.se_acc).^2;
            costmat(:,:,:,j,k,7) = ((squeeze(data.b1_RT(:,:,k,:)) - coeff_rt_sub_or.beta_OV)/coeff_rt_sub_or.se_OV).^2;
            costmat(:,:,:,j,k,8) = ((squeeze(data.b2_RT(:,:,k,:)) - coeff_rt_sub_or.beta_VD)/coeff_rt_sub_or.se_VD).^2;
            costmat(:,:,:,j,k,9) = ((squeeze(data.b1_Acc(:,:,k,:)) - coeff_acc_sub_or.beta_OV)/coeff_acc_sub_or.se_OV).^2;
            costmat(:,:,:,j,k,10) = ((squeeze(data.b2_Acc(:,:,k,:)) - coeff_acc_sub_or.beta_VD)/coeff_acc_sub_or.se_VD).^2;
            costmat(:,:,:,j,k,11) = ((bestndt + squeeze(data.I_RT(:,:,k,:)) - intercept_sub_or.I_rt)/intercept_sub_or.se_rt).^2;
            costmat(:,:,:,j,k,12) = ((squeeze(data.I_Acc(:,:,k,:)) - intercept_sub_or.I_acc)/intercept_sub_or.se_acc).^2;
            bestndtmat(:,:,:,j,k) = bestndt;
        end    
    end
    cost = sum(costmat,6);
    if min(cost(:)) < mincost
        % mincost(j) = min(cost(:));
        [v,loc] = min(cost(:));
        mincost = v;
        [ii,jj,kk,ww,pp] = ind2sub(size(cost),loc);
        for m = 1:2
            decays(m) = data.decay(ii,jj,j,kk);
            mis(m) = data.mutualInhib(ii,jj,j,kk);
            gains(m) = data.gain(ii,jj,j,kk);
            orders(m) = data.order;
            nDTs(m) = data.nDT;
            nondecsT(m) = bestndtmat(ii,jj,kk,ww,pp);
        end
        
        threshs(1) = data.thresh(ii,jj,ww,kk);
        threshs(2) = data.thresh(ii,jj,pp,kk);
        file = files(i).name;
    end
end


%%

besttheta = [gains(1) 0 decays(1) decays(1)*mis(1) threshs(1) nondecsT(1) 0.01 0.5;
             gains(1) 0 decays(1) decays(1)*mis(1) threshs(2) nondecsT(1) 0.01 0.5];

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
%%
table_outcome = array2table(outcome');

table_outcome.Properties.VariableNames(1:6) = {'I_RT','b1_RT','b2_RT','I_Acc','b1_Acc','b2_Acc'};

table_outcome.condition = {'XOR';'OR'};

writetable(table_outcome,'../result/prediction_task_12_thresh.csv');