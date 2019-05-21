function dyna = dynet_sim(n,Fs,duration,order,sparsity,fdom,nstates,ntrials)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Simulation framework for tv-MVAR generated surrogate time series
%                                       D. Pascucci, University of Fribourg
% Last update: 18.05.2019 - D. Pascucci
%--------------------------------------------------------------------------
% INPUT:
% n        (number of nodes)      Fs      (sampling freq)    
% duration (trial lenght, sec)    order   (model order)    
% sparsity (proportion)           fdom    (univariate freq, Hz)
% nstates  (number of states)     ntrials (number of realizations)
%--------------------------------------------------------------------------
% OUTPUT (dynet, structure with fields)
% dyne.net     sub-structure with the main properties of surrogate data
% dyna.SC      binary structural connections matrix
% dyna.DC      binary asymmetric functional (directed) connectivity matrix
% dyna.AR      time varying AR coefficients (n x n x order x time)
% dyna.Y       time-series (trials x n x time)
% dyna.CT      across-trials correlation matrix
% dyna.NS      noise covariance matrix (Identity for the moment)
% dyna.scaling scaling factor for off-diagonal elements (for stability)
% dyna.regimes temporal segments defying network states

%--------------------------------------------------------------------------
% TODO:
% connectivity summary in a table
% add linear mixing?
% add snr?
% noise covariance?
% FURTHER IMPROVEMENTS:
% 3D lattice of realistic distances for SC
% small-world properties
% time-varying covariance and mean (ERP)
% projection onto the scalp with a realistic head model

%--------------------------------------------------------------------------
% defaults                                                                 
default('n',5);        default('Fs',200); 
default('duration',2); default('order',15);   
default('sparsity',.5);default('fdom',[8 16]);
default('nstates',2);  default('ntrials',200);
global srate p time 
head       = {'state','rec','send','band','mag','time','osc'};
%--------------------------------------------------------------------------
% constants
SClinks    = 0.8;           % Markov et al., 2012
ascale     = 0.05:0.01:0.3; % 0.2 max better
srate      = Fs;
dt         = 1/Fs;
time       = 0:dt:(duration-dt);
samples    = numel(time);
p          = order;
lag_t      = dt:dt:(p*dt);
min_state  = unique(dsearchn(time',duration*.2)); % in frames
%--------------------------------------------------------------------------
% structural links (binary mask)
I          = eye(n);
UT         = triu(reshape(randsample([0 1],n^2,'true',...
                       [1-SClinks SClinks]),[n n]));
MK         = (UT+UT') - diag(diag(UT+UT'));
SC         = MK + I;
%--------------------------------------------------------------------------
% directed interactions (binary mask)
DC         = zeros(size(SC));
DC(randsample(find(MK),fix((1-sparsity)*numel(find(MK))),'false')) = 1;
DC         = DC + I;
%--------------------------------------------------------------------------
% ar process (univariate)
AR         = zeros(n,n,p,samples);
ampli      = randsample(ascale,n,'true');
OSC        = damposc(repmat(fdom,[n 1]),lag_t+(dt/2),0.05,ampli,(p*dt)*.5);% @Joan: fix the lag_t+(dt/2) issue
for i = 1:n
    AR(i,i,:,:) = repmat(OSC(i,:)',[1 1 numel(time)]);
end
%--------------------------------------------------------------------------
% ar process (interactions)
cON        = find(MK.*DC);
start_at   = unique(dsearchn(time',...
                    randsample(duration*.1:dt:duration*.2,1)));
state_ons  = sort(randsample(start_at:min_state:(samples-min_state), ...
                                                            nstates));
state_end  = state_ons+randsample(min_state:samples,nstates);
state_dur  = min(state_end,diff([state_ons samples]));
% determine states and check stability
regimes{nstates} = [];
scalef           = 1;
summary_conn     = table();
for k = 1:nstates
    regimes{k}   = state_ons(k):state_ons(k)+state_dur(k);
    ok           = 0;
    summary      = NaN(numel(cON),8+p);
    while ok==0
        % generate off-diag AR and check stability
        tmpAR      = AR(:,:,:,regimes{k}(1));
        for ij = 1:numel(cON)
            freq   = sort(randsample(5:fix(Fs/2)-5,2));
            ampli  = randsample(ascale,1)*scalef;            
            osc    = damposc(freq(1):freq(2),lag_t+(dt/2),0.05,...
                                                          ampli,(p*dt)*.5); 
            [i,j]  = ind2sub([n n],cON(ij));
            summary(ij,:) = [k i j freq ampli ...
                       time(regimes{k}(1)) time(regimes{k}(end)) ...
                                    osc'];
            tmpAR(i,j,:)  = osc;
        end
        blockA     = [tmpAR(:,:); eye((p-1)*n) zeros((p-1)*n,n)];
        if any(abs(eig(blockA))>.95)
            scalef = scalef*.95;
        else
            ok     = 1;
        end
    end
    tmp            = table(summary(:,1),summary(:,2),summary(:,3),...
                           summary(:,4:5),summary(:,6),summary(:,7:8),...
                           summary(:,9:end),'VariableNames',head);                       
    summary_conn   = vertcat(summary_conn,tmp);
    % add stable matrices to the dynamical system
    AR(:,:,:,regimes{k})  = repmat(tmpAR,[1 1 1 numel(regimes{k})]);
end
% AR          = movingmean(AR,0.2/dt,4);
%--------------------------------------------------------------------------
% data (add nuisance segment at the beginning)
X           = zeros(ntrials,n,numel(time)+start_at);
ARplus      = cat(4,AR(:,:,:,1:start_at),AR);
CT          = min(abs(gallery('randcorr',ntrials))*3,1);
dgI         = shuffling(find(eye(ntrials)==0));
CT(dgI(1:fix(numel(dgI)*.1))) = -CT(dgI(1:fix(numel(dgI)*.1)));
% C           = diag(max(abs(randsample(-3:0.5:3,n,'true')),.5));
for k_p = 1:p
    X(:,:,k_p)    = CT*randn(ntrials,n,1);
end
for k = (p+1):numel(time)+start_at
    for l = 1:p
        X(:,:,k)  =  X(:,:,k) + (ARplus(:,:,l,k)*X(:,:,k-l)')' + ...
                                    CT*randn(ntrials,n,1);% * C;           % how about the covariance?
    end
end
X(:,:,1:start_at) = []; % remove nuisance data
AR                = AR(:,:,:,1:samples); % ensure size
%==========================================================================
dyna.net     = struct('n',n,'srate',Fs,'p',p,'time',time,   ...            % store freq of interaction?
                      'sparsity',sparsity,'regimes',nstates,...
                      'trials',ntrials);
dyna.SC      = SC;
dyna.DC      = DC;
dyna.AR      = AR;
dyna.Y       = X;
dyna.CT      = CT;
dyna.R       = eye(n);% C;
dyna.scaling = scalef;
dyna.regimes = regimes;
dyna.summary = summary_conn;

% clc
% disp('Simulated dynamic functional network')
% disp(dyna.net)
