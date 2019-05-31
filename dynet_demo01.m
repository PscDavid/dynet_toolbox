free
addpath(genpath('E:\Git\tvConnectivity\Dynet_tool\dynet_toolbox'))
%--------------------------------------------------------------------------
% Simulation
%--------------------------------------------------------------------------
% simulate surrogate data with default network settings 
% (type 'help dynet_sim' for customizing) 
sim        = dynet_sim(5,200,2,5,.5,3,400,-1);
% review the main properties of the simulated network
review(sim)
% compute and visualize the ground truth PDC obtained directly from the 
% true tvMVAR data-generating process (here squared-PDC, column norm)
gt_PDC     = dynet_ar2pdc(sim,sim.srate,sim.frange','sPDC',[],[],1);
dynet_connplot(gt_PDC,sim.time,sim.frange',[],[],[],sim.DC,1)
% note: in the diagonal the parametric MVAR-derived PSD, scaled to the
%       range of off-diagonal PDC values for graphical purpose;
%       cells framed in red are open (existing) functional connections

%--------------------------------------------------------------------------
% Adaptive filtering
%--------------------------------------------------------------------------
% general Linear Kalman Filter (canonical c=0.01)
KF         = dynet_SSM_KF(sim.Y,sim.popt,0.01);
kf_PDC     = dynet_ar2pdc(KF,sim.srate,sim.frange,'sPDC',[],[],1);
dynet_connplot(kf_PDC,sim.time,sim.frange,[],[],[],sim.DC,1)
% Sparse Adaptive Least-squares Kalman
% lambda     = dynet_ISS_GCV(sim.Y,sim.popt,[],[],1);  % determine lambda GCV
SA         = dynet_SSM_SALKff(sim.Y,sim.popt,0.99);
sa_PDC     = dynet_ar2pdc(SA,sim.srate,sim.frange,'sPDC',[],[],1);
dynet_connplot(sa_PDC,sim.time,sim.frange,[],[],[],sim.DC,1)

%--------------------------------------------------------------------------
% Compare
%--------------------------------------------------------------------------
% compare the two algorithms with ROC analysis
figure
auc_kf     = roc_auc(gt_PDC,kf_PDC,20,1);hold on
auc_sa     = roc_auc(gt_PDC,sa_PDC,20,1);
disp(['AUC results: KF = ' num2str(auc_kf) ' SALK = ' num2str(auc_sa)])
