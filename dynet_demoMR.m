clc; clear all; close all;

addpath('E:\Git\tvConnectivity\Dynet_tool\dynet_toolbox\connectivity')
addpath('E:\Git\tvConnectivity\Dynet_tool\dynet_toolbox\simulation')
addpath('E:\Git\tvConnectivity\Dynet_tool\dynet_toolbox\statespace')
addpath('E:\Git\tvConnectivity\Dynet_tool\dynet_toolbox\statespace\optimization')
addpath('E:\Git\tvConnectivity\Dynet_tool\dynet_toolbox\utilities')
addpath('E:\Git\tvConnectivity\Dynet_tool\dynet_toolbox\utilities\external')
addpath('E:\Git\tvConnectivity\Dynet_tool\dynet_toolbox\utilities\external\cbrewer\cbrewer')

%%
global srate p time
%--------------------------------------------------------------------------
% Simulation
%--------------------------------------------------------------------------
% simulate surrogate data with default network settings 
% (type 'help dynet_sim' for customizing) 
sim        = dynet_sim();
% review the main properties of the simulated network
review(sim)
% compute and visualize the ground truth PDC obtained directly from the 
% true tvMVAR data-generating process (here squared-PDC, column norm)
freqs      = (1:srate/2)'; % frequency range of interest 
gt_PDC     = dynet_ar2pdc(sim,srate,freqs,'sPDC',[],[],1);
dynet_connplot(gt_PDC,time,freqs,[],[],[],sim.DC,1)
% note: in the diagonal is the parametric MVAR-derived PSD, scaled to the
%       range of the off-diagonal PDC for graphic purpose;
%       cells framed in red are open (existing) funcitonal connections

%--------------------------------------------------------------------------
% Adaptive filtering
%--------------------------------------------------------------------------
% general Linear Kalman Filter (canonical c=0.01)
KF         = dynet_SSM_KF(sim.Y,p,0.01);
kf_PDC     = dynet_ar2pdc(KF,srate,freqs,'sPDC',[],[],1);
dynet_connplot(kf_PDC,time,freqs,[],[],[],sim.DC,1)
% Sparse Adaptive Least-squares Kalman
lambda     = dynet_ISS_GCV(sim.Y,p,[],[],1);  % determine lambda GCV
SA         = dynet_SSM_SALK(sim.Y,p,lambda);
sa_PDC     = dynet_ar2pdc(SA,srate,freqs,'sPDC',[],[],1);
dynet_connplot(sa_PDC,time,freqs,[],[],[],sim.DC,1)

%--------------------------------------------------------------------------
% Compare
%--------------------------------------------------------------------------
% compare the two algorithms with ROC analysis
figure
auc_kf     = roc_auc(gt_PDC,kf_PDC,20,1);hold on
auc_sa     = roc_auc(gt_PDC,sa_PDC,20,1);
legend([' KF=' num2str(auc_kf)],'',['SALK=' num2str(auc_sa)])
disp(['AUC results: KF = ' num2str(auc_kf) ' SALK = ' num2str(auc_sa)])
