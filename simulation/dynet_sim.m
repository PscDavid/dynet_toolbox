function dyna = dynet_sim(n,Fs,duration,order,sparsity,...
                          nstates,ntrials,snr_db,lmix)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Simulation framework for tv-MVAR generated surrogate time series
%                                       
% Last update: 23.08.2019
%--------------------------------------------------------------------------
% INPUTs:
% - n:       number of nodes           - Fs:     sampling frequency
% - duration:trial lenght in s         - order:  model order
% - sparsity:proportion                - nstates:number of states
% - ntrials: number of realizations    - snr_db: signal-to-noise ratio db
% - lmix:    linear mixing in percentage
%--------------------------------------------------------------------------
% OUTPUT dyna, structure with fields:
% - n:       number of nodes             - srate:  sampling frequency Fs
% - delay:   model order                 - popt:   optimal model order
% - time:    time vector                 - frange: frequency vector
% - sparsity:matrix sparsity             - nstates:number of states
% - trials:  number of realizations      - SC:     structural adjacency matrix
% - DC:      functional adjacency matrix - AR:     MVAR model coefficients
% - R:       Innovation covariance       - Y:      Simulated sample data
% - E:       Innovation matrix           - CT:     Between-trials correlation
% - scaling: Scale factor of off-diag AR - regimes:Stable matrices
% - DM,LMx:  Distance matrix             - summary:Table of connections
%--------------------------------------------------------------------------
% INVOKED FUNCTIONs
% default.m; shuffling.m; addnoise.m
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% - defaults
default('n',5);         default('Fs',200);
default('duration',2);  default('order',fix(0.025/(1/Fs)));
default('sparsity',.5); default('nstates',3);
default('ntrials',200); default('snr_db',NaN);
default('lmix',0);
head       = {'state','rec','send','mag','time','osc','lagop'};

% - constants
SClinks    = 0.8;                         % Markov et al., 2012
ascale     = 0.1:0.01:0.5;
dt         = 1/Fs;
time       = 0:dt:(duration-dt);
samples    = numel(time);
p          = order;
min_state  = unique(dsearchn(time',.15)); % minimum duration in frames

% - structural links (binary mask)
I          = eye(n);
UT         = triu(reshape(randsample([0 1],n^2,'true',...
    [1-SClinks SClinks]),[n n]));
MK         = (UT+UT') - diag(diag(UT+UT'));
SC         = MK + I;

% - directed interactions (binary mask)
DC         = zeros(size(SC));
DC(randsample(find(MK),fix((1-sparsity)*numel(find(MK))),'false')) = 1;
DC         = DC + I;

% - AR process (univariate)
AR         = zeros(n,n,p+1,samples);
for i = 1:n
    c1            = randsample(ascale,1);
    c2            = (max(ascale)-c1)*.95;
    AR(i,i,1:2,:) = repmat([c1 c2],[1 1 numel(time)]); % low-pass
end

% - AR process (interactions)
cON        = find(MK.*DC);
[bf,~]     = buffer(1:numel(time),min_state);
start_at   = sort(randsample(3:size(bf,2),nstates));
state_ons  = bf(1,start_at);
state_end  = [bf(1,start_at(2:end)) numel(time)];                          % no empty states in between?
state_dur  = state_end-state_ons;
% determine states
regimes{nstates} = [];
% starting (no scaling)
scalef           = 1;
summary_conn     = table();
for k = 1:nstates
    regimes{k}   = state_ons(k):state_ons(k)+state_dur(k);
    ok           = 0;
    summary      = NaN(numel(cON),9);
    while ok==0
        % generate off-diag AR and check stability
        tmpAR      = AR(:,:,:,regimes{k}(1));
        for ij = 1:numel(cON)
            ij_p   = randsample(1:p,1);
            ampl1  = randsample(ascale*.5,1);                              % scaling off diag?
            ampl2  = (max(ascale*.5)-ampl1)*0.95;
            osc    = sign(randn(1,2)).*[ampl1 ampl2].*scalef;
            [i,j]  = ind2sub([n n],cON(ij));
            summary(ij,:) = [k i j ampl1 ...
                time(regimes{k}(1)) time(regimes{k}(end)) ...
                osc ij_p];
            tmpAR(i,j,ij_p:ij_p+1)  = osc;
        end
        pp         = p+1;
        % stability check
        blockA     = [tmpAR(:,:); eye((pp-1)*n) zeros((pp-1)*n,n)];
        if any(abs(eig(blockA))>.95)
            scalef = scalef*.95;
        else
            ok     = 1;
        end
    end
    tmp            = table(summary(:,1),summary(:,2),summary(:,3),...
        summary(:,4),summary(:,5:6),...
        summary(:,7:end-1),summary(:,end),...
        'VariableNames',head);
    summary_conn   = vertcat(summary_conn,tmp);
    % add stable matrices to the dynamical system
    AR(:,:,:,regimes{k})  = repmat(tmpAR,[1 1 1 numel(regimes{k})]);
end
% AR          = movingmean(AR,0.2/dt,4);

% - Data (add nuisance segment at the beginning)
nuisance    = min(state_ons(1)*2,.5/dt);
X           = zeros(ntrials,n,numel(time)+nuisance);
ARplus      = cat(4,AR(:,:,:,1:nuisance),AR);
% simulate between-trials correlation (correlated generative noise)
% CT          = triu(reshape(normrnd(0.3,0.1,ntrials^2,1),[ntrials ntrials]));
% CT          = CT+CT' - diag(diag(CT+CT')) + eye(ntrials);
% CT          = min(CT,1);
% CT=eye(ntrials);
CT          = min(abs(gallery('randcorr',ntrials))*3,1); %3
dgI         = shuffling(find(eye(ntrials)==0));
CT(dgI(1:fix(numel(dgI)*.1))) = -CT(dgI(1:fix(numel(dgI)*.1)));

% - Generate time-series
for k_p = 1:size(AR,3)                     % the actual p or popt
    X(:,:,k_p)    = CT*randn(ntrials,n,1);
end
E           = X;
for k = (size(AR,3)+1):numel(time)+nuisance
    innovation    = CT*randn(ntrials,n,1); % across trials correlation
    E(:,:,k)      = innovation;
    X(:,:,k)      = X(:,:,k)+innovation;
    for l = 1:size(AR,3)
        X(:,:,k)  = X(:,:,k) + (ARplus(:,:,l,k)*X(:,:,k-l)')';
    end
end
X(:,:,1:nuisance) = [];                    % remove nuisance data
Y                 = X;                     % observed signal
E(:,:,1:nuisance) = [];                    % innovation
AR                = AR(:,:,:,1:samples);   % ensure size, AR coeffs
tmp               = permute(E,[2 3 1]);
R                 = (tmp(:,:)*tmp(:,:)');  % innovation noise covariance

% - SNR
if ~isnan(snr_db)
    Y              = addnoise(Y,snr_db,'w');
end

% - Linear mixing
if lmix>0
    x             = randsample(1:1:150,n,'false')'; % 2d lattice 15x15 cm
    y             = randsample(1:1:150,n,'false')';
    mixf          = normpdf(0:150,0,lmix);
    xy            = [x y];
    [j,i]         = meshgrid(1:n,1:n);
    distance      = xy(i,:) - xy(j,:);
    DM            = zeros(n);
    DM(:)         = sqrt(sum(distance.^2,2));       % distance matrix (cm)
    LMx           = normpdf(DM,0,lmix);
    LMx           = LMx./max(mixf);
    for tr = 1:ntrials
        Y(tr,:,:)     = LMx*squeeze(Y(tr,:,:));
    end
else
    DM            = [];
    LMx           = [];
end

% - Output structure
dyna              = struct('n',n,'srate',Fs,'delay',p,             ...
    'popt',size(AR,3),'time',time,          ...
    'frange',(1:fix(Fs/2))',                ...
    'sparsity',sparsity,'nstates',nstates,  ...
    'trials',ntrials,'SC',SC,'DC',DC,       ...
    'AR',AR,'Y',Y,'E',E,'CT',CT,'R',R,      ...
    'scaling',scalef,'regimes',{regimes},   ...
    'DM',DM,'LMx',LMx,                      ...
    'summary',summary_conn);

end
