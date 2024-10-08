function [summary] = summaryPredictionSeparate(thetas,deadline,order)


M = readtable('../../data/trial/simulation_value_12.csv');
inputs = sort(M{:,:},2,'descend');


% Is = [];
% 
% for v1 = 2:10
%     for v2 = 1:(v1-1)
%         for v3 = 1:v2
%             for v4 = 1:v3
%                 Is = [Is;[v1,v2,v3,v4]];
%             end
%         end
%     end
% end
% inputs = sort(Is,2,'descend');
% 
nTrials = size(inputs,1);

meanValue = mean(inputs,2);
valueDiff = inputs(:,1) - mean(inputs(:,[2 3 4]),2);

%% ===== Set up hyper-parameters =====

%%

% r = [0 0.5];

summary = [];

for j = 1:2
    theta = thetas(j,:);
    
    options.doRelu = true;
    options.nDT    = deadline(j);
    options.dt     = 0.05;
    options.maxRT  = options.nDT * options.dt;     % maximum RT
    options.nRep   = 100;     % number of repititions of process
    options.bound  = (1 - linspace(0,1,options.nDT+1).^order(j));    % boundary (1 by (nDT + 1))

    RTlist = [];
    Hitlist = [];
    meanValuelist = [];
    valueDifflist = [];

    for m = 1:nTrials

        [RTs,Choices] = SSM_collapse_gen(theta,inputs(m,:),options);
        RTlist = [RTlist;RTs];
        Hitlist = [Hitlist;Choices==1];
        meanValuelist = [meanValuelist;repmat(meanValue(m),options.nRep,1)];
        valueDifflist = [valueDifflist;repmat(valueDiff(m),options.nRep,1)];
    end
    
    % figure;hist(RTlist);

    rtfit = fitlm(zscore([meanValuelist valueDifflist]),...
        RTlist);
    accfit = fitglm(zscore([meanValuelist valueDifflist]),...
        Hitlist,'Distribution','binomial');

    I_RT = rtfit.Coefficients.Estimate(1);
    b1_RT = rtfit.Coefficients.Estimate(2);
    b2_RT = rtfit.Coefficients.Estimate(3);
    I_Acc = exp(accfit.Coefficients.Estimate(1))/(1+exp(accfit.Coefficients.Estimate(1)));
    b1_Acc = accfit.Coefficients.Estimate(2);
    b2_Acc = accfit.Coefficients.Estimate(3);
    
    summary(1,j) = I_RT;
    summary(2,j) = b1_RT;
    summary(3,j) = b2_RT;
    summary(4,j) = I_Acc;
    summary(5,j) = b1_Acc;
    summary(6,j) = b2_Acc;
end


end

