clear all; close all;

files = dir('../outcome/*.mat');

%%
tic;
datalist = {};

for i = 1:length(files)
    datalist{i} = load(['../outcome/' files(i).name]);
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
bestfiles = {};
coeff_rt_sub_or = coeff_rt(strcmp(coeff_rt.study,study) & strcmp(coeff_rt.catChoiceCondition,'OR'),:);
coeff_acc_sub_or = coeff_acc(strcmp(coeff_acc.study,study) & strcmp(coeff_acc.catChoiceCondition,'OR'),:);
intercept_sub_or = intercept(strcmp(intercept.study,study) & strcmp(intercept.catChoiceCondition,'OR'),:);
coeff_rt_sub_xor = coeff_rt(strcmp(coeff_rt.study,study) & strcmp(coeff_rt.catChoiceCondition,'XOR'),:);
coeff_acc_sub_xor = coeff_acc(strcmp(coeff_acc.study,study) & strcmp(coeff_acc.catChoiceCondition,'XOR'),:);
intercept_sub_xor = intercept(strcmp(intercept.study,study) & strcmp(intercept.catChoiceCondition,'XOR'),:);

%% Fit XOR

mintotalcost = 100000;

for t = 0.05:0.05:1
    for m = 1:11
        for k = 1:11
            mincost = 100000;

            totalcost = 0;

            for i = 1:length(datalist)
                data = datalist{i};% load(['../outcome/' files(i).name]);

                costmat = zeros(11,11,6);

                % bestndtmat = intercept_sub_xor.I_rt - squeeze(data.I_RT(:,m,:,k));
                % bestndtmat(bestndtmat<0.5) = 0.5;
                costmat(:,:,1) = ((squeeze(data.b1_RT(:,m,:,k)) - coeff_rt_sub_xor.beta_OV)/coeff_rt_sub_xor.se_OV).^2;
                costmat(:,:,2) = ((squeeze(data.b2_RT(:,m,:,k)) - coeff_rt_sub_xor.beta_VD)/coeff_rt_sub_xor.se_VD).^2;
                costmat(:,:,3) = ((squeeze(data.b1_Acc(:,m,:,k)) - coeff_acc_sub_xor.beta_OV)/coeff_acc_sub_xor.se_OV).^2;
                costmat(:,:,4) = ((squeeze(data.b2_Acc(:,m,:,k)) - coeff_acc_sub_xor.beta_VD)/coeff_acc_sub_xor.se_VD).^2;
                costmat(:,:,5) = ((t + squeeze(data.I_RT(:,m,:,k)) - intercept_sub_xor.I_rt)/intercept_sub_xor.se_rt).^2;
                costmat(:,:,6) = ((squeeze(data.I_Acc(:,m,:,k)) - intercept_sub_xor.I_acc)/intercept_sub_xor.se_acc).^2;
            %             costmat(:,:,:,j,k,7) = ((squeeze(data.b1_RT(k,:,:,:)) - coeff_rt_sub_or.beta_OV)/coeff_rt_sub_or.se_OV).^2;
            %             costmat(:,:,:,j,k,8) = ((squeeze(data.b2_RT(k,:,:,:)) - coeff_rt_sub_or.beta_VD)/coeff_rt_sub_or.se_VD).^2;
            %             costmat(:,:,:,j,k,9) = ((squeeze(data.b1_Acc(k,:,:,:)) - coeff_acc_sub_or.beta_OV)/coeff_acc_sub_or.se_OV).^2;
            %             costmat(:,:,:,j,k,10) = ((squeeze(data.b2_Acc(k,:,:,:)) - coeff_acc_sub_or.beta_VD)/coeff_acc_sub_or.se_VD).^2;
            %             costmat(:,:,:,j,k,11) = 0 * ((bestndt + squeeze(data.I_RT(k,:,:,:)) - intercept_sub_or.I_rt)/intercept_sub_or.se_rt).^2;
            %             costmat(:,:,:,j,k,12) = 0 * ((squeeze(data.I_Acc(k,:,:,:)) - intercept_sub_or.I_acc)/intercept_sub_or.se_acc).^2;

                %cost = max(costmat,[],3);
                cost = sum(costmat,3);
                if min(cost(:)) < mincost
                    % mincost(j) = min(cost(:));
                    [v,loc] = min(cost(:));
                    mincost = v;
                    [ii,jj] = ind2sub(size(cost),loc);
                    decays(1) = data.decay(ii,m,jj,k);
                    threshs(1) = data.thresh(ii,m,jj,k);
                    gains(1) = data.gain(ii,m,jj,k);
                    orders(1) = data.order;
                    nDTs(1) = data.nDT;
                    nondecsT(1) = t;
                    mis(1) = data.mutualInhib(ii,m,jj,k);
                    bestfiles{1} = files(i).name;
                end
            end

            totalcost = mincost;

            mincost = 100000;
            for i = 1:length(datalist)
                data = datalist{i};% load(['../outcome/' files(i).name]);

                costmat = zeros(11,11,6);

                % bestndtmat = intercept_sub_or.I_rt - squeeze(data.I_RT(:,m,:,k));
                % bestndtmat(bestndtmat<0.5) = 0.5;
                
                
                
                
                costmat(:,:,1) = ((squeeze(data.b1_RT(:,m,:,k)) - coeff_rt_sub_or.beta_OV)/coeff_rt_sub_or.se_OV).^2;
                costmat(:,:,2) = ((squeeze(data.b2_RT(:,m,:,k)) - coeff_rt_sub_or.beta_VD)/coeff_rt_sub_or.se_VD).^2;
                costmat(:,:,3) = ((squeeze(data.b1_Acc(:,m,:,k)) - coeff_acc_sub_or.beta_OV)/coeff_acc_sub_or.se_OV).^2;
                costmat(:,:,4) = ((squeeze(data.b2_Acc(:,m,:,k)) - coeff_acc_sub_or.beta_VD)/coeff_acc_sub_or.se_VD).^2;
                costmat(:,:,5) = ((t + squeeze(data.I_RT(:,m,:,k)) - intercept_sub_or.I_rt)/intercept_sub_or.se_rt).^2;
                costmat(:,:,6) = ((squeeze(data.I_Acc(:,m,:,k)) - intercept_sub_or.I_acc)/intercept_sub_or.se_acc).^2;
            %     costmat(:,:,:,j,k,7) = ((squeeze(data.b1_RT(k,:,:,:)) - coeff_rt_sub_or.beta_OV)/coeff_rt_sub_or.se_OV).^2;
            %     costmat(:,:,:,j,k,8) = ((squeeze(data.b2_RT(k,:,:,:)) - coeff_rt_sub_or.beta_VD)/coeff_rt_sub_or.se_VD).^2;
            %     costmat(:,:,:,j,k,9) = ((squeeze(data.b1_Acc(k,:,:,:)) - coeff_acc_sub_or.beta_OV)/coeff_acc_sub_or.se_OV).^2;
            %     costmat(:,:,:,j,k,10) = ((squeeze(data.b2_Acc(k,:,:,:)) - coeff_acc_sub_or.beta_VD)/coeff_acc_sub_or.se_VD).^2;
            %     costmat(:,:,:,j,k,11) = 0 * ((bestndt + squeeze(data.I_RT(k,:,:,:)) - intercept_sub_or.I_rt)/intercept_sub_or.se_rt).^2;
            %     costmat(:,:,:,j,k,12) = 0 * ((squeeze(data.I_Acc(k,:,:,:)) - intercept_sub_or.I_acc)/intercept_sub_or.se_acc).^2;
                cost = sum(costmat,3);
                if min(cost(:)) < mincost
                    [v,loc] = min(cost(:));
                    mincost = v;
                    [ii,jj,kk] = ind2sub(size(cost),loc);
                    decays(2) = data.decay(ii,m,jj,k);
                    threshs(2) = data.thresh(ii,m,jj,k);
                    gains(2) = data.gain(ii,m,jj,k);
                    orders(2) = data.order;
                    nDTs(2) = data.nDT;
                    nondecsT(2) = t;
                    mis(2) = data.mutualInhib(ii,m,jj,k);
                    bestfiles{2} = files(i).name;
                end
            end

            totalcost = totalcost + mincost;

            if totalcost < mintotalcost

                mintotalcost = totalcost;
                bestdecays = decays;
                bestthreshs = threshs;
                bestgains = gains;
                bestorders = orders;
                bestnDTs = nDTs;
                bestnondecsT = nondecsT;
                bestmis = mis;
                bestbestfiles = bestfiles;

            end

        end
    end
end
%%

besttheta = [bestgains(1) 0 bestdecays(1) bestdecays(1)*bestmis(1) bestthreshs(1) bestnondecsT(1) 0.01 0.5;
             bestgains(2) 0 bestdecays(2) bestdecays(2)*bestmis(2) bestthreshs(2) bestnondecsT(2) 0.01 0.5];

outcome = summaryPredictionSeparate(besttheta,bestnDTs,bestorders);

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

writetable(table_outcome,'../result/prediction_task_12_full.csv');

save("../result/params_full.mat","besttheta","bestnDTs","bestorders");