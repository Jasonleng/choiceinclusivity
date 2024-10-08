function [RT,Choice,finalState] = SSM_collapse_gen(theta,input,options)

% ===== Parameter of interest for fitting =====

forward     = theta(1);          % gain of input
inhibit     = theta(2);          % feedforward inhibition of input
decay       = theta(3);          % decay
competition = theta(4);          % competition
thresh      = theta(5);          % decision threshold
t0          = theta(6);          % non-decision time
lapse       = theta(7);          % lapse rate
sd          = theta(8);          % diffusion noise

% ===== simulation parameters  =====

doRelu = options.doRelu;    % force positive accumulators
nDT    = options.nDT;       % max number of timesteps
dt     = options.dt;        % time step
maxRT  = options.maxRT;     % maximum RT
nRep   = options.nRep;      % number of repititions of process
bound  = options.bound;     % boundary (1 by (nDT + 1))

% Set up collapsing boundary

bound  = thresh * bound;

% ===== Size of simulation =====

nChoices = length(input);          % how many accumulators?

% ===== Setup noise sampling parameters =====

sddt = sd*sqrt(dt);

% ===== Setup dynamic matrix =====

a_lca = eye(nChoices) - (competition * ones(nChoices) + (decay - competition)*eye(nChoices))*dt;

a_ffi = -inhibit * ones(nChoices) + (forward + inhibit) * eye(nChoices);

% ===== Setup input column =====

b = a_ffi * input' * dt;

% Initialize matrix for storing simulation

x = zeros(nChoices,nRep); % Current status stored in x
all_x = nan(nDT+1,nChoices,nRep); % History stored in all_x
converge = max(x) > thresh; % Status of respone stored in converge

% ===== Run simulation =====

for tt = 1:nDT+1
    
    all_x(tt,:,:) = x;
    
    converge = converge | (max(x)>bound(tt));
    
    if sum(converge) == nRep
        break;
    end
    
    x = a_lca * x + b + randn(size(x))*sddt;
    
    % Nonlinear function to avoid negative x value
    
    if doRelu
        x = max(x,0);
    end
    
end

% ===== Extract RT and Choice =====

[RT, Choice] = deal(nan(nRep, 1));
finalState = nan(nRep,nChoices);

for rr = 1:nRep
    try
        RT(rr)      = find(max(all_x(:,:,rr),[],2)>bound', 1);
        finalState(rr,:)  = all_x(RT(rr),:,rr);
        [~,Choice(rr)]  = max(all_x(RT(rr),:,rr));
    catch
        RT(rr) = nan;
        finalState(rr,:)  = all_x(nDT+1,:,rr);
        Choice(rr) = 0;
    end
end

% ===== Include lapse trials =====

lapse_idx = rand(size(RT))<lapse;

lapse_Choice = randi(nChoices,[sum(lapse_idx),1]);

Choice(lapse_idx) = lapse_Choice;

RT = (RT-1)*dt; % Transfer RT from steps to seconds

lapse_RT = rand([sum(lapse_idx),1]) * maxRT;

RT(lapse_idx) = lapse_RT;

RT = RT + t0; % Include ndt;

end

