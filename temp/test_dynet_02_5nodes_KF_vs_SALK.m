free
mainfold   = 'C:\Users\David\Google Drive\Work\09_Code\DYN_Tool\dynet_toolbox';
addpath(genpath(mainfold))
%--------------------------------------------------------------------------
% General Linear Kalman Filter vs. Sparse Adaptive Least-squares Filter
% ROC analysis
%--------------------------------------------------------------------------
% - Basic 5 nodes network, varying c
c_grid     = logspace(-4,0,10);
nsim       = 20;
n          = 50;
% - preallocate
auc_sa     = NaN(nsim,1);      c_sa = auc_sa;
auc_kf     = NaN(nsim,numel(c_grid));
% - loop over simulations
for b = 1:nsim
    sim         = dynet_sim(n);
    gt_PDC      = dynet_ar2pdc(sim,sim.srate,sim.frange,'sPDC',0,2);
    % estimate SALK model
    SA          = dynet_SSM_SALKff(sim.Y,sim.popt,0.99);
    sa_PDC      = dynet_ar2pdc(SA,sim.srate,sim.frange,'sPDC',0,2);
    % SALK auc
    auc_sa(b,:) = roc_auc(gt_PDC,sa_PDC,20);
    c_sa(b,:)   = mean(SA.c);
    % loop over c for general Linear Kalman Filter 
    for k = 1:numel(c_grid)
        c              = c_grid(k);
        KF             = dynet_SSM_KF(sim.Y,sim.popt,c);
        kf_PDC         = dynet_ar2pdc(KF,sim.srate,sim.frange,'sPDC',0,2);
        auc_kf(b,k)    = roc_auc(gt_PDC,kf_PDC,20);
    end
    disp(['simulation ' num2str(b) ' of ' num2str(nsim)])
end
% - save
cd([mainfold '\test\results']);
save(['results_dynet_02_' num2str(n) 'nodes_KF_vs_SALK.mat'], ...
      'c_grid', 'c_sa', 'auc_sa', 'auc_kf', 'nsim')

%% plot results
% Kalman filtering over adaptation constant
figure('position',[488  195  389 566])
hp  = ciplot(c_grid,auc_kf,[87 108 67]./255,0);
set(gca,'xscale','log','xlim',c_grid([1 end]),'ylim',[.5 1]);
xlabel('\itc');ylabel('Area under the ROC');
tl  = title('Kalman filter');tl.FontWeight = 'normal';
figformat(0,0,.1)
% Sparse Adaptive Least-squares Kalman with self-tuning c
figure('position',[488  195  195 566])
hp2 = ciplot(mean(c_sa),auc_sa,[220 0 0]./255,0);
ciplot(mean(c_sa)+[-.1 0 .1],repmat(auc_sa,[1 3]),[220 0 0]./255);
xlim(mean(c_sa)+[-.15 .15]);ylim([.5 1])
set(gca,'xtick',mean(c_sa),'xticklabel',num2str(mean(c_sa),'%.2f'))
xlabel('average \itc');
tl  = title('SALK filter');tl.FontWeight = 'normal';
figformat(0,0,.1)
