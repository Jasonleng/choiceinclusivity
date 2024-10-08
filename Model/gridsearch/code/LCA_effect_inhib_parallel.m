function LCA_effect_inhib_parallel(jid,tid)
%% ===== Value of products =====
addpath('../../models');
M = readtable('../../data/trial/simulation_value_12.csv');
inputs = sort(M{:,1:4},2,'descend');
nTrials = size(inputs,1);

%% ===== Read from inputs and assign jobs =====

if tid < 121

nDT   = fix(tid/11)*10+60;
order = (rem(tid,11))*0.1+1;

%% ===== Set up hyper-parameters =====

options.doRelu = true;
options.nDT    = nDT;
options.dt     = 0.05;
options.maxRT  = options.nDT * options.dt;     % maximum RT
options.nRep   = 100;     % number of repititions of process
options.bound  = (1 - linspace(0,1,options.nDT+1).^order);    % boundary (1 by (nDT + 1))

%% ===== Get some IVs =====

meanValue = mean(inputs,2);
valueDiff = inputs(:,1) - mean(inputs(:,[2 3 4]),2);

%% ===== Set up levels of parameters =====

% Intensity of decay

dMin = 0;
dMax = 5;
dNlevel = 11;
dLevels = linspace(dMin,dMax,dNlevel);

% Ratio between competition and decay

mMin = 0;
mMax = 2;
mNlevel = 41;
mLevels = linspace(mMin,mMax,mNlevel);

% threshold

aMin = 0.5;
aMax = 2.5;
aNlevel = 11;
aLevels = linspace(aMin,aMax,aNlevel);

% intensity of signal

cMin = 0.05;
cMax = 0.55;
cNlevel = 11;
cLevels = linspace(cMin,cMax,cNlevel);

%% Simulate the effect of decay and inhibition

% Set up holder for beta values

I_RT = nan(mNlevel,dNlevel,aNlevel,cNlevel);
I_Acc = nan(mNlevel,dNlevel,aNlevel,cNlevel);
b1_RT = nan(mNlevel,dNlevel,aNlevel,cNlevel);
b1_Acc = nan(mNlevel,dNlevel,aNlevel,cNlevel);
b2_RT = nan(mNlevel,dNlevel,aNlevel,cNlevel);
b2_Acc = nan(mNlevel,dNlevel,aNlevel,cNlevel);
decay = nan(mNlevel,dNlevel,aNlevel,cNlevel);
mutualInhib = nan(mNlevel,dNlevel,aNlevel,cNlevel);
thresh = nan(mNlevel,dNlevel,aNlevel,cNlevel);
gain = nan(mNlevel,dNlevel,aNlevel,cNlevel);

tic;

for k1 = 1:dNlevel
    for k2 = 1:aNlevel
        for k3 = 1:cNlevel
            for k4 = 1:mNlevel
                disp(['=====> Starting fitting # ', num2str(k1), '-',num2str(k2),'-',num2str(k3),'-',num2str(k4)]);

                theta_true = [cLevels(k3) 0 dLevels(k1) dLevels(k1)*mLevels(k4) aLevels(k2) 0 0.01 0.5];
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

                I_RT(k4,k1,k2,k3) = rtfit.Coefficients.Estimate(1);
                b1_RT(k4,k1,k2,k3) = rtfit.Coefficients.Estimate(2);
                b2_RT(k4,k1,k2,k3) = rtfit.Coefficients.Estimate(3);
                I_Acc(k4,k1,k2,k3) = exp(accfit.Coefficients.Estimate(1))/(1+exp(accfit.Coefficients.Estimate(1)));
                b1_Acc(k4,k1,k2,k3) = accfit.Coefficients.Estimate(2);
                b2_Acc(k4,k1,k2,k3) = accfit.Coefficients.Estimate(3);
                decay(k4,k1,k2,k3) = dLevels(k1);
                mutualInhib(k4,k1,k2,k3) = mLevels(k4);
                thresh(k4,k1,k2,k3) = aLevels(k2);
                gain(k4,k1,k2,k3) = cLevels(k3);
            end
        end
    end
end


toc;
save(sprintf('../outcome/LCA_%d_task_%d.mat',jid,tid),'nDT','order','decay','mutualInhib','thresh','gain','I_Acc','I_RT','b1_Acc','b1_RT','b2_Acc','b2_RT');

end

end
