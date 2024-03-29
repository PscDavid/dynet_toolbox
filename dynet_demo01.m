%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Brief tutorial of dynet_toolbox
%
% Last update 22.10.2019
%--------------------------------------------------------------------------
% INVOKED FUNCTIONs
% free.m; dynet_sim.m; review.m; dynet_ar2pdc.m; dynet_connplot.m; 
% dynet_SSM_KF.m; dynet_SSM_STOK.m; roc_auc.m
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%  Initializing

% I1) Add the path containing the folder "dynet_toolbox" 

mypath     = ['C:\Users\chare\Google Drive\Work\09_Code\DYN_Tool\dynet_toolbox'];
p          = genpath(mypath);
addpath(p)

% I2) Clear everything

free

%% Simulation

% S1) Simulate surrogate data with default network settings 
% (type 'help dynet_sim' for customizing) 

sim        = dynet_sim(10,200,2,5,.6,3,400,5,0);

% S2) Review the main properties of the simulated network

review(sim)

% S3) Compute and visualize the ground truth PDC obtained directly from the 
% true tvMVAR data-generating process (here squared-PDC, column norm)

gt_PDC     = dynet_ar2pdc(sim,sim.srate,sim.frange','sPDC',0,[],1);
dynet_connplot(gt_PDC,sim.time,sim.frange',[],[],[],sim.FC,1)

% note: in the diagonal the parametric MVAR-derived PSD, scaled to the
%       range of off-diagonal PDC values for graphical purpose;
%       cells framed in red are open (existing) functional connections

%% Adaptive Filtering

% AF1) Apply the general Linear Kalman Filter (canonical c=0.02) & compute
% sPDC

KF         = dynet_SSM_KF(sim.Y,sim.popt,0.02);
kf_PDC     = dynet_ar2pdc(KF,sim.srate,sim.frange,'sPDC',[],[],1);
dynet_connplot(kf_PDC,sim.time,sim.frange,[],[],[],sim.FC,1)
sgtitle('gKF') % only with >MATLAB R2019a

% AF2) Apply the STOK filter & compute sPDC

SK         = dynet_SSM_STOK(sim.Y,sim.popt);
sk_PDC     = dynet_ar2pdc(SK,sim.srate,sim.frange,'sPDC',[],[],1);
dynet_connplot(sk_PDC,sim.time,sim.frange,[],[],[],sim.FC,1)
sgtitle('STOK') % only with >MATLAB R2019a

%% Comparison

% C1) Compare the two algorithms with ROC analysis

figure()
auc_kf     = roc_auc(gt_PDC,kf_PDC,20,1);
hold on
auc_sk     = roc_auc(gt_PDC,sk_PDC,20,1);
legend('gKF','','STOK')
disp(['AUC results: KF = ' num2str(auc_kf) ' STOK = ' num2str(auc_sk)])
