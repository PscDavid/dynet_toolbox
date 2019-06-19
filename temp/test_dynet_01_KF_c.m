free
mainfold   = 'C:\Users\David\Google Drive\Work\09_Code\DYN_Tool\dynet_toolbox';
addpath(genpath(mainfold))
%--------------------------------------------------------------------------
% General Linear Kalman Filter
%--------------------------------------------------------------------------
% - small vs large network, varying c
c_grid     = logspace(-3,0,10);
n_grid     = [10 100];
nsim       = 1;
% - preallocate
aucKF      = NaN(numel(c_grid),numel(n_grid),nsim);
crKF       = aucKF;
% - loop over simulations
for b = 1:nsim
    sim         = dynet_sim(10);
    gt_PDC      = dynet_ar2pdc(sim,sim.srate,sim.frange,'sPDC',0,2);
    nonanel     = ~isnan(gt_PDC);
    for n = 1:numel(n_grid) % add noisy channels
        noise     = normrnd(mean(sim.Y(:)),std(sim.Y(:)),...
                     sim.trials,n_grid(n)-sim.n,numel(sim.time));
        Y         = cat(2,sim.Y,noise);
        for k = 1:numel(c_grid) % loop over c
            c            = c_grid(k);
            KF           = dynet_SSM_KF(Y,sim.popt,c);
            kf_PDC       = dynet_ar2pdc(KF,sim.srate,sim.frange,'sPDC',0,2);
            kf_PDC       = kf_PDC(1:sim.n,1:sim.n,:,:);
            aucKF(k,n,b) = roc_auc(gt_PDC(nonanel),kf_PDC(nonanel),20);
            crKF(k,n,b)  = corr(gt_PDC(nonanel),kf_PDC(nonanel));

            disp(['sim ' num2str(b) ' n ' num2str(n_grid(n)) ' c ' ...
                num2str(c)])
        end
    end
end
% plot(c_grid,aucKF);set(gca,'xscale','log')
cd(mainfold);cd('test\results');
save results_dynet_01_KF_c.mat c_grid n_grid aucKF crKF nsim

%% ========================================================================
% % variance
% c_grid     = logspace(-3,0,10);
% sim        = dynet_sim(10);
% % preallocate
% varFreqC   = NaN(numel(c_grid),1);
% varC       = NaN(sim.n*(sim.n-1)*numel(sim.frange)*numel(sim.time),...
%                  numel(c_grid));
% for k = 1:numel(c_grid)
%     c             = c_grid(k);
%     KF            = dynet_SSM_KF(sim.Y,sim.popt,c);
%     kf_PDC        = dynet_ar2pdc(KF,sim.srate,sim.frange,'sPDC',0,2);
%     varFreqC(k,:) = nansum(nansum(nansum(nanvar(kf_PDC,[],3))));
%     varC(:,k)     = kf_PDC(~isnan(kf_PDC));
% end
% delta_c    = sum(diff(varC,[],2).^2);
% cd(mainfold);cd('test\results');
% save(ls('*01_KF_c.mat'),'varFreqC','varC','delta_c','-append')
% 
% subplot(1,2,1);
% plot(c_grid,varFreqC,'linewidth',3);set(gca,'xscale','log');
% xlabel('\itc');ylabel('spectral variability');xlim(c_grid([1 end]));
% figformat(0,0,0.1);
% 
% subplot(1,2,2);
% plot(c_grid(2:end),delta_c,'linewidth',3);set(gca,'xscale','log');
% xlabel('\itc');ylabel('\deltaPDC');xlim(c_grid([2 end]));
% figformat(0,0,0.1);

% %% ========================================================================
% % additional (p order variance)
% p_grid     = 2:2:30;
% sim        = dynet_sim(10);
% % preallocate
% varFreqP   = NaN(numel(p_grid),1);
% varP       = NaN(sim.n*(sim.n-1)*numel(sim.frange)*numel(sim.time),...
%                  numel(p_grid));
% for k = 1:numel(p_grid)
%     p             = p_grid(k);
%     KF            = dynet_SSM_KF(sim.Y,p,0.04);
%     kf_PDC        = dynet_ar2pdc(KF,sim.srate,sim.frange,'sPDC',0,2);
%     varFreqP(k,:) = nansum(nansum(nansum(nanvar(kf_PDC,[],3))));
%     varP(:,k)     = kf_PDC(~isnan(kf_PDC));
% end
% delta_p    = sum(diff(varP,[],2).^2);
% cd(mainfold);cd('test\results');
% % save(ls('*01_KF_c.mat'),'varFreqP','varP','delta_p','-append')
% 
% subplot(1,2,1);
% plot(p_grid,varFreqP,'linewidth',3);
% xlabel('\itp');ylabel('spectral variability');xlim(p_grid([1 end]));
% figformat(0,0,0.1);
% 
% subplot(1,2,2);
% plot(p_grid(2:end),delta_p,'linewidth',3);
% xlabel('\itp');ylabel('\deltaPDC');xlim(p_grid([2 end]));
% figformat(0,0,0.1);
