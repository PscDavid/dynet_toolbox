%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%                    General Linear Kalman Filter (KF)
%                                  VS.
%               Self-tuning Optimized Kalman filter (STOK)
%                                       D. Pascucci, University of Fribourg
% Last update: 17.06.2019
%--------------------------------------------------------------------------
% Metrics:
% 1) Area Under the ROC (AUC) over c (KF) 
% 2) PDC correlation over p
% 3) AUC, over SNR
% 4) AUC, number of trials * nodes
% 5) AUC, over linear mixing
% 6) Speed [regression] over nodes
%==========================================================================
free
mfold  = 'C:\Users\David\Google Drive\Work\09_Code\DYN_Tool\dynet_toolbox';
addpath(genpath(mfold))

%% 1) Area Under the ROC (AUC) over c (KF) 
%--------------------------------------------------------------------------
% main variables
nsim       = 30;
clist      = logspace(-4,0,10);
%--------------------------------------------------------------------------
% network variables
nodes      = 10;
%--------------------------------------------------------------------------
% preallocate
AUC_kf     = NaN(nsim,numel(clist));
AUC_sk     = NaN(nsim,1);
%--------------------------------------------------------------------------
% main loop
for i = 1:nsim
    % simulate ground truth
    sim         = dynet_sim(nodes);
    PDC_gt      = dynet_ar2pdc(sim,sim.srate,sim.frange,'sPDC',[],2);
    % STOK
    SK          = dynet_SSM_STOK(sim.Y,sim.popt,.99);
    PDC_sk      = dynet_ar2pdc(SK,sim.srate,sim.frange,'sPDC',[],2);
    AUC_sk(i,:) = roc_auc(PDC_gt,PDC_sk,20,0);
    % nested loop over clist
    for k = 1:numel(clist)
        c           = clist(k);
        % KF
        KF          = dynet_SSM_KF(sim.Y,sim.popt,c);
        PDC_kf      = dynet_ar2pdc(KF,sim.srate,sim.frange,'sPDC',[],2);
        AUC_kf(i,k) = roc_auc(PDC_gt,PDC_kf,20,0);
    end
    prcdone(i,nsim,'step 1) AUC over c.',10)
end
cd(mfold)
cd('test\results\simulated_data')
save KF_vs_STOK_1_c_constant.mat clist nsim AUC_kf AUC_sk

%% 2) PDC correlation over p
%--------------------------------------------------------------------------
% main variables
nsim       = 30;
plist      = 2:15;
coptKF     = 0.0167; % choose from step 1
%--------------------------------------------------------------------------
% network variables
nodes      = 10;
fs         = 200;
tlength    = 2;
samples    = tlength/(1/fs);
p          = 5;    % plus 2 for the addition of on-zero AR(2)
%--------------------------------------------------------------------------
% preallocate
AUC_kf     = NaN(nsim,numel(plist));
AUC_sk     = AUC_kf;
CORR_KF    = NaN(numel(plist),numel(plist),nsim);
CORR_SK    = CORR_KF;
%--------------------------------------------------------------------------
% main loop
for i = 1:nsim
    % simulate ground truth
    sim         = dynet_sim(nodes,fs,p);
    PDC_gt      = dynet_ar2pdc(sim,sim.srate,sim.frange,'sPDC',[],2);
    % nested loop over plist
    PDC_VC_kf   = NaN(nnz(~isnan(PDC_gt)),numel(plist));
    PDC_VC_sk   = PDC_VC_kf;
    for k = 1:numel(plist)
        pk          = plist(k);
        % STOK
        SK          = dynet_SSM_STOK(sim.Y,pk,.99);
        PDC_sk      = dynet_ar2pdc(SK,sim.srate,sim.frange,'sPDC',[],2);
        PDC_VC_sk(:,k) = PDC_sk(~isnan(PDC_gt));
        AUC_sk(i,k) = roc_auc(PDC_gt,PDC_sk,20,0);
        % KF
        KF          = dynet_SSM_KF(sim.Y,pk,coptKF);
        PDC_kf      = dynet_ar2pdc(KF,sim.srate,sim.frange,'sPDC',[],2);
        PDC_VC_kf(:,k) = PDC_kf(~isnan(PDC_gt));
        AUC_kf(i,k) = roc_auc(PDC_gt,PDC_kf,20,0);
    end
    CORR_KF(:,:,i)  = corr(PDC_VC_kf);
    CORR_SK(:,:,i)  = corr(PDC_VC_sk);
    prcdone(i,nsim,'step 2) PDC correlation over p.',10)
end
cd(mfold)
cd('test\results\simulated_data')
save KF_vs_STOK_2_porder.mat plist nsim AUC_kf AUC_sk CORR_KF CORR_SK

%% 3) AUC, over SNR
%--------------------------------------------------------------------------
% main variables
nsim       = 30;
noise      = [0.1,1,3,5,10];
coptKF     = 0.0167; % choose from step 1
%--------------------------------------------------------------------------
% network variables
nodes      = 10;
%--------------------------------------------------------------------------
% preallocate
AUC_kf     = NaN(nsim,numel(noise));
AUC_sk     = AUC_kf;
%--------------------------------------------------------------------------
% main loop
for n = 1:numel(noise)
    for i = 1:nsim
        sim         = dynet_sim(nodes,[],[],[],[],[],[],noise(n));
        PDC_gt      = dynet_ar2pdc(sim,sim.srate,sim.frange,'sPDC',[],2);
        % KF
        KF          = dynet_SSM_KF(sim.Y,sim.popt,coptKF); 
        PDC_kf      = dynet_ar2pdc(KF,sim.srate,sim.frange,'sPDC',[],2);
        AUC_kf(i,n) = roc_auc(PDC_gt,PDC_kf,20,0);
        % STOK
        SK          = dynet_SSM_STOK(sim.Y,sim.popt,.99);
        PDC_sk      = dynet_ar2pdc(SK,sim.srate,sim.frange,'sPDC',[],2);
        AUC_sk(i,n) = roc_auc(PDC_gt,PDC_sk,20,0);       
    end
    prcdone(n,numel(noise),'step 3) AUC, over SNR.',10)
end
cd(mfold)
cd('test\results\simulated_data')
save KF_vs_STOK_3_snr.mat noise nsim AUC_kf AUC_sk

%% 4) AUC, number of trials * nodes
%--------------------------------------------------------------------------
% main variables
nsim       = 30;
ntrials    = [100 200 300 400 500];
ndim       = [5 10 20 40];
coptKF     = 0.0167; % choose from step 1
%--------------------------------------------------------------------------
% preallocate
AUC_kf     = NaN(nsim,numel(ndim),numel(ntrials));
AUC_sk     = AUC_kf;
%--------------------------------------------------------------------------
% main loop
for tr = 1:numel(ntrials)
    for n = 1:numel(ndim)
        for i = 1:nsim
            sim       = dynet_sim(ndim(n),[],[],[],[],[],ntrials(tr),[]);
            PDC_gt    = dynet_ar2pdc(sim,sim.srate,sim.frange,'sPDC',[],2);
            % kalman
            KF        = dynet_SSM_KF(sim.Y,sim.popt,coptKF); 
            PDC_kf    = dynet_ar2pdc(KF,sim.srate,sim.frange,'sPDC',[],2);
            AUC_kf(i,n,tr) = roc_auc(PDC_gt,PDC_kf,20,0);
            % regularized
            SK        = dynet_SSM_STOK(sim.Y,sim.popt,.99);
            PDC_sk    = dynet_ar2pdc(SK,sim.srate,sim.frange,'sPDC',[],2);
            AUC_sk(i,n,tr) = roc_auc(PDC_gt,PDC_sk,20,0); 
        end
    end
    prcdone(tr,numel(ntrials),'step 4) AUC, number of trials * nodes.',10)
end
cd(mfold)
cd('test\results\simulated_data')
save KF_vs_STOK_4_trialsXnodes.mat ntrials ndim nsim AUC_kf AUC_sk

%% 5) AUC, over linear mixing
%--------------------------------------------------------------------------
% main variables
nsim       = 30;
sp         = linspace(0.1,0.5,5);
coptKF     = 0.0167; % choose from step 1
%--------------------------------------------------------------------------
% network variables
nodes      = 10;
%--------------------------------------------------------------------------
% preallocate
AUC_kf     = NaN(nsim,numel(sp));
AUC_sk     = AUC_kf;
%--------------------------------------------------------------------------
% main loop
for n = 1:numel(sp)
    for i = 1:nsim
        sim         = dynet_sim(nodes,[],[],[],[],[],[],[],sp(n));
        PDC_gt      = dynet_ar2pdc(sim,sim.srate,sim.frange,'sPDC',[],2);
        % kalman
        KF          = dynet_SSM_KF(sim.Y,sim.popt,coptKF); 
        PDC_kf      = dynet_ar2pdc(KF,sim.srate,sim.frange,'sPDC',[],2);
        AUC_kf(i,n) = roc_auc(PDC_gt,PDC_kf,20,0);
        % regularized
        SK          = dynet_SSM_STOK(sim.Y,sim.popt,.99);
        PDC_sk      = dynet_ar2pdc(SK,sim.srate,sim.frange,'sPDC',[],2);
        AUC_sk(i,n) = roc_auc(PDC_gt,PDC_sk,20,0); 
    end
    prcdone(n,numel(sp),'step 5) AUC, over linear mixing.',10)
end
cd(mfold)
cd('test\results\simulated_data')
save KF_vs_STOK_5_mixing.mat sp nsim AUC_kf AUC_sk

%% 6) Speed [regression] over nodes
%--------------------------------------------------------------------------
% main variables
nsim       = 30;
ndim       = [10 40];
%--------------------------------------------------------------------------
% preallocate
R2KF       = NaN(nsim,numel(ndim),100);
R2SK       = R2KF;
for i = 1:nsim
    for n = 1:numel(ndim)
        sim         = dynet_sim(ndim(n),[],[],[],[],[],200);
        PDC_gt      = dynet_ar2pdc(sim,sim.srate,sim.frange,'sPDC',[],2);
        % KF
        KF          = dynet_SSM_KF(sim.Y,sim.popt,0.017); 
        PDC_kf      = dynet_ar2pdc(KF,sim.srate,sim.frange,'sPDC',[],2);
        % STOK
        SK          = dynet_SSM_STOK(sim.Y,sim.popt,.99);
        PDC_sk      = dynet_ar2pdc(SK,sim.srate,sim.frange,'sPDC',[],2);
    %----------------------------------------------------------------------
    % nested loops
    r2_kf    = NaN(numel(sim.frange),numel(sim.regimes)-1);
    r2_sk    = r2_kf; delta  = r2_kf;
    for r = 1:numel(sim.regimes)-1
        xf          = cat(4,PDC_gt(:,:,:,sim.regimes{r}(1)  +(-2:2)),...
                            PDC_gt(:,:,:,sim.regimes{r}(end)+(-2:2)));
        yf_kf       = cat(4,PDC_kf(:,:,:,sim.regimes{r}(1)  +(-2:2)),...
                            PDC_kf(:,:,:,sim.regimes{r}(end)+(-2:2)));
        yf_sk       = cat(4,PDC_sk(:,:,:,sim.regimes{r}(1)  +(-2:2)),...
                            PDC_sk(:,:,:,sim.regimes{r}(end)+(-2:2)));
        for fr = 1:numel(sim.frange)
            % KF
            x       = squeeze(xf(:,:,fr,:));    
            y       = squeeze(yf_kf(:,:,fr,:)); 
            x       = zscore(x(~isnan(y))); y       = zscore(y(~isnan(y)));
            r2_kf(fr,r) = 1- sum(square( y-(y'/x')*x) ) ./ sum(square( y ) );
            % STOK
            y       = squeeze(yf_sk(:,:,fr,:));
            y       = zscore(y(~isnan(y)));
            r2_sk(fr,r) = 1- sum(square( y-(y'/x')*x) ) ./ sum(square( y ) );
        end
    end
    % store
    R2KF(i,n,:) = mean(r2_kf,2);
    R2SK(i,n,:) = mean(r2_sk,2);
    end
    prcdone(i,nsim,'6) Speed [regression] over nodes.',10)
end
cd(mfold)
cd('test\results\simulated_data')
save KF_vs_STOK_6_speed_over_freq.mat ndim nsim R2KF R2SK

%%
%%=========================================================================
% plots
free
mfold  = 'C:\Users\David\Google Drive\Work\09_Code\DYN_Tool\dynet_toolbox';
addpath(genpath(mfold))
% 1) Area Under the ROC (AUC) over c (KF) 
load(ls('*1*'))
% cols   = [0.3373    0.7059    0.8745; 0.8353    0.3686    0];
cols   = [0.6    0.6    0.6; 0.8353    0.3686    0];
cp     = ciplot(clist,AUC_kf,cols(1,:));hold on
set(gca,'xscale','log');
[m,ci] = estimates(0.05, AUC_sk, 0);
pl     = plot(0.14,m,'.','color',cols(2,:),'markersize',40);hold on
eb     = errorbar(0.14,m,m-ci(1,:),ci(2,:)-m); eb.Color = [0 0 0];
xlabel('\it c');ylabel('AUC');ylim([.6 1])
formatlegend([cp pl],{'KF','STOK'})
figformat
tl     = title('KF vs. STOK (c)');tl.FontWeight = 'normal';
%--------------------------------------------------------------------------
% 2) PDC correlation over p
clearvars -except mfold cols
load(ls('*2*'))
figure
subplot(1,3,1);% KF
imagesc(plist,plist,mean(CORR_KF,3));axis xy;
tl     = title('KF(p corr)');tl.FontWeight = 'normal';
caxis([0.5 1]);xlabel('p');ylabel('p')
figformat;cb = colorbar; cb.LineWidth = 1.5;
subplot(1,3,2);% STOK
imagesc(plist,plist,mean(CORR_SK,3));axis xy;
tl     = title('STOK(p corr)');tl.FontWeight = 'normal';
caxis([0.5 1]);xlabel('p');ylabel('p');
figformat;cb = colorbar; cb.LineWidth = 1.5;
subplot(1,3,3);% DIFFERENCE
imagesc(plist,plist,mean(CORR_SK-CORR_KF,3));axis xy;
tl     = title('\delta (STOK-KF)');tl.FontWeight = 'normal';
caxis([-.2 .2]);xlabel('p');ylabel('p');
figformat;cb = colorbar; cb.LineWidth = 1.5;
%--------------------------------------------------------------------------
% 3) AUC, over SNR
clearvars -except mfold cols
load(ls('*3*'))
figure
cp(1)  = ciplot(noise,AUC_kf,cols(1,:));hold on
cp(2)  = ciplot(noise,AUC_sk,cols(2,:));hold on
set(gca,'xscale','log');
xlabel('snr(db)');ylabel('AUC');ylim([.6 1])
figformat
formatlegend(cp,{'KF','STOK'},5,20)
tl     = title('KF vs. STOK (snr)');tl.FontWeight = 'normal';
%--------------------------------------------------------------------------
% 4) AUC, number of trials * nodes
clearvars -except mfold cols
load(ls('*4*'))
figure
imagesc(squeeze(mean(AUC_sk-AUC_kf)));
set(gca,'xtick',1:5,'xticklabels',ntrials,'ytick',1:4,'yticklabels',ndim)
caxis([-.15 .15]);cb = colorbar; cb.LineWidth = 1.5;
ylabel('nodes');xlabel('trials');
figformat
%--------------------------------------------------------------------------
% 5) AUC, over linear mixing
clearvars -except mfold cols
load(ls('*5*'))
figure
cp(1)  = ciplot(sp,AUC_kf,cols(1,:));hold on
cp(2)  = ciplot(sp,AUC_sk,cols(2,:));hold on
% set(gca,'xscale','log');
xlabel('linear mixing');ylabel('AUC');ylim([.6 1])
figformat
formatlegend(cp,{'KF','STOK'},5,20)
tl     = title('KF vs. STOK (mixing)');tl.FontWeight = 'normal';
%--------------------------------------------------------------------------
% 6) Speed [regression] over nodes
clearvars -except mfold cols
load(ls('*6*'))
figure
% n = 10
subplot(1,2,1)
tmp    = squeeze(R2KF(:,1,:));
cp(1)  = ciplot(1:100,tmp,cols(1,:));hold on
tmp    = squeeze(R2SK(:,1,:));
cp(2)  = ciplot(1:100,tmp,cols(2,:));hold on
xlabel('frequency(Hz)');
ylim([0.1 0.8]);ylabel('R^{2}');
tl     = title('n: 10; o:200');tl.FontWeight = 'normal';
figformat
formatlegend(cp,{'KF','STOK'},5,20)
% n = 40
subplot(1,2,2)
tmp    = squeeze(R2KF(:,2,:));
cp(1)  = ciplot(1:100,tmp,cols(1,:));hold on
tmp    = squeeze(R2SK(:,2,:));
cp(2)  = ciplot(1:100,tmp,cols(2,:));hold on
xlabel('frequency(Hz)');
ylim([0.1 0.8]);ylabel('R^{2}')
tl     = title('n: 40; o:200');tl.FontWeight = 'normal';
figformat

