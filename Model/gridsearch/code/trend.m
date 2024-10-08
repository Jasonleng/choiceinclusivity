clear all;close all;


%% ===== Value of products =====
addpath('../../models');
M = readtable('../../data/trial/simulation_value_12.csv');
inputs = sort(M{:,1:4},2,'descend');
nTrials = size(inputs,1);


%%

params = load('../result/params_inhib.mat');

%% ===== Set up hyper-parameters =====

options.doRelu = true;
options.nDT    = params.nDTs(2);
options.dt     = 0.05;
options.maxRT  = options.nDT * options.dt;     % maximum RT
options.nRep   = 200;     % number of repititions of process
options.bound  = (1 - linspace(0,1,options.nDT+1).^params.orders(2));    % boundary (1 by (nDT + 1))

%%

meanValue = mean(inputs,2);
valueDiff = inputs(:,1) - mean(inputs(:,[2 3 4]),2);

%%

mMin = 0;
mMax = 1;
mNlevel = 11;
mLevels = linspace(mMin,mMax,mNlevel);


I_RT = nan(mNlevel,1);
b1_RT = nan(mNlevel,1);
b2_RT = nan(mNlevel,1);
I_Acc = nan(mNlevel,1);
b1_Acc = nan(mNlevel,1);
b2_Acc = nan(mNlevel,1);
mutualInhib = nan(mNlevel,1);

tic;

for k4 = 1:mNlevel
    % disp(['=====> Starting fitting # ', num2str(k1), '-',num2str(k2),'-',num2str(k3),'-',num2str(k4)]);

    theta_true = params.besttheta(1,:);
    theta_true(4) = theta_true(3) * mLevels(k4);
    RTlist = nan(nTrials * options.nRep,1);
    Hitlist = nan(nTrials * options.nRep,1);
    meanValuelist = nan(nTrials * options.nRep,1);
    valueDifflist = nan(nTrials * options.nRep,1);

    for j = 1:nTrials
        [RTs,Choices] = SSM_collapse_gen(theta_true,inputs(j,:),options);
        RTlist((1 + (j - 1) * options.nRep):(j * options.nRep),1) = RTs;
        Hitlist((1 + (j - 1) * options.nRep):(j * options.nRep),1) = (Choices==1);
        meanValuelist((1 + (j - 1) * options.nRep):(j * options.nRep),1) = repmat(meanValue(j),options.nRep,1);
        valueDifflist((1 + (j - 1) * options.nRep):(j * options.nRep),1) = repmat(valueDiff(j),options.nRep,1);
    end

    rtfit = fitlm(zscore([meanValuelist valueDifflist]),...
        RTlist);
    accfit = fitglm(zscore([meanValuelist valueDifflist]),...
        Hitlist,'Distribution','binomial');

    I_RT(k4) = rtfit.Coefficients.Estimate(1);
    b1_RT(k4) = rtfit.Coefficients.Estimate(2);
    b2_RT(k4) = rtfit.Coefficients.Estimate(3);
    I_Acc(k4) = exp(accfit.Coefficients.Estimate(1))/(1+exp(accfit.Coefficients.Estimate(1)));
    b1_Acc(k4) = accfit.Coefficients.Estimate(2);
    b2_Acc(k4) = accfit.Coefficients.Estimate(3);
    mutualInhib(k4) = mLevels(k4);
end

%%

outcome = [mutualInhib I_RT b1_RT b2_RT I_Acc b1_Acc b2_Acc];

table_outcome = array2table(outcome);

table_outcome.Properties.VariableNames(1:7) = {'m','I_RT','b1_RT','b2_RT','I_Acc','b1_Acc','b2_Acc'};

writetable(table_outcome,'../result/trend_m.csv');

%%
aMin = 0.4;
aMax = 1.7;
aNlevel = 20;
aLevels = linspace(aMin,aMax,aNlevel);


I_RT = nan(aNlevel,1);
b1_RT = nan(aNlevel,1);
b2_RT = nan(aNlevel,1);
I_Acc = nan(aNlevel,1);
b1_Acc = nan(aNlevel,1);
b2_Acc = nan(aNlevel,1);
threshold = nan(aNlevel,1);

tic;

for k4 = 1:aNlevel
    % disp(['=====> Starting fitting # ', num2str(k1), '-',num2str(k2),'-',num2str(k3),'-',num2str(k4)]);

    theta_true = params.besttheta(1,:);
    theta_true(5) = aLevels(k4);
    RTlist = nan(nTrials * options.nRep,1);
    Hitlist = nan(nTrials * options.nRep,1);
    meanValuelist = nan(nTrials * options.nRep,1);
    valueDifflist = nan(nTrials * options.nRep,1);

    for j = 1:nTrials
        [RTs,Choices] = SSM_collapse_gen(theta_true,inputs(j,:),options);
        RTlist((1 + (j - 1) * options.nRep):(j * options.nRep),1) = RTs;
        Hitlist((1 + (j - 1) * options.nRep):(j * options.nRep),1) = (Choices==1);
        meanValuelist((1 + (j - 1) * options.nRep):(j * options.nRep),1) = repmat(meanValue(j),options.nRep,1);
        valueDifflist((1 + (j - 1) * options.nRep):(j * options.nRep),1) = repmat(valueDiff(j),options.nRep,1);
    end

    rtfit = fitlm(zscore([meanValuelist valueDifflist]),...
        RTlist);
    accfit = fitglm(zscore([meanValuelist valueDifflist]),...
        Hitlist,'Distribution','binomial');

    I_RT(k4) = rtfit.Coefficients.Estimate(1);
    b1_RT(k4) = rtfit.Coefficients.Estimate(2);
    b2_RT(k4) = rtfit.Coefficients.Estimate(3);
    I_Acc(k4) = exp(accfit.Coefficients.Estimate(1))/(1+exp(accfit.Coefficients.Estimate(1)));
    b1_Acc(k4) = accfit.Coefficients.Estimate(2);
    b2_Acc(k4) = accfit.Coefficients.Estimate(3);
    threshold(k4) = aLevels(k4);
end

%%

outcome = [threshold I_RT b1_RT b2_RT I_Acc b1_Acc b2_Acc];

table_outcome = array2table(outcome);

table_outcome.Properties.VariableNames(1:7) = {'a','I_RT','b1_RT','b2_RT','I_Acc','b1_Acc','b2_Acc'};

writetable(table_outcome,'../result/trend_a.csv');

%%
dMin = 60;
dMax = 160;
dNlevel = 11;
dLevels = linspace(dMin,dMax,dNlevel);


I_RT = nan(dNlevel,1);
b1_RT = nan(dNlevel,1);
b2_RT = nan(dNlevel,1);
I_Acc = nan(dNlevel,1);
b1_Acc = nan(dNlevel,1);
b2_Acc = nan(dNlevel,1);
ddl = nan(dNlevel,1);

tic;

for k4 = 1:dNlevel
    % disp(['=====> Starting fitting # ', num2str(k1), '-',num2str(k2),'-',num2str(k3),'-',num2str(k4)]);

    theta_true     = params.besttheta(1,:);
    options.nDT    = dLevels(k4);
    options.maxRT  = options.nDT * options.dt;     % maximum RT
    options.bound  = (1 - linspace(0,1,options.nDT+1).^params.orders(2));    % boundary (1 by (nDT + 1))    
    
    RTlist = nan(nTrials * options.nRep,1);
    Hitlist = nan(nTrials * options.nRep,1);
    meanValuelist = nan(nTrials * options.nRep,1);
    valueDifflist = nan(nTrials * options.nRep,1);

    for j = 1:nTrials
        [RTs,Choices] = SSM_collapse_gen(theta_true,inputs(j,:),options);
        RTlist((1 + (j - 1) * options.nRep):(j * options.nRep),1) = RTs;
        Hitlist((1 + (j - 1) * options.nRep):(j * options.nRep),1) = (Choices==1);
        meanValuelist((1 + (j - 1) * options.nRep):(j * options.nRep),1) = repmat(meanValue(j),options.nRep,1);
        valueDifflist((1 + (j - 1) * options.nRep):(j * options.nRep),1) = repmat(valueDiff(j),options.nRep,1);
    end

    rtfit = fitlm(zscore([meanValuelist valueDifflist]),...
        RTlist);
    accfit = fitglm(zscore([meanValuelist valueDifflist]),...
        Hitlist,'Distribution','binomial');

    I_RT(k4) = rtfit.Coefficients.Estimate(1);
    b1_RT(k4) = rtfit.Coefficients.Estimate(2);
    b2_RT(k4) = rtfit.Coefficients.Estimate(3);
    I_Acc(k4) = exp(accfit.Coefficients.Estimate(1))/(1+exp(accfit.Coefficients.Estimate(1)));
    b1_Acc(k4) = accfit.Coefficients.Estimate(2);
    b2_Acc(k4) = accfit.Coefficients.Estimate(3);
    ddl(k4) = dLevels(k4);
end

%%

outcome = [ddl I_RT b1_RT b2_RT I_Acc b1_Acc b2_Acc];

table_outcome = array2table(outcome);

table_outcome.Properties.VariableNames(1:7) = {'theta','I_RT','b1_RT','b2_RT','I_Acc','b1_Acc','b2_Acc'};

writetable(table_outcome,'../result/trend_theta.csv');