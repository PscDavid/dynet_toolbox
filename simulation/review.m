function review(sim)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Figure displaying: 
% suplot(221) structural adjacency matrix
% suplot(222) functional adjacency matrix
% subplot(223) surrogate time-series in the time domain
% subplot(224) power spectral density of surrogate time-series
%
% Last update 22.08.2019
%--------------------------------------------------------------------------
% INVOKED FUNCTIONs
% ShadePlotForEmpahsis.m; figformat.m
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

figure('units','normalized','position',[0.25 0.1 .5 .8]);
subplot(2,2,1)
imagesc(sim.SC);
colormap(gca,[0 0 0; [0.8500 0.3250 0.0980]]);
% plot(1:sim.net.n,1:sim.net.n,'k-');
set(gca,'xtick',linspace(1,sim.n,sim.n)-.5,'xticklabel',[],...
    'xgrid','on','xcolor','w',...
    'ytick',linspace(1,sim.n,sim.n)-.5,'ytickLabel',[],...
    'ygrid','on','ycolor','w',...
    'gridLineStyle','-','linewidth',2,'gridalpha',1)
xlabel('nodes'); ylabel('nodes');
figformat;
tl = title('Structural connections');tl.FontWeight = 'normal';
subplot(2,2,2)
imagesc(sim.FC);
colormap(gca,[0 0 0; [0.8500 0.3250 0.0980]]);
% plot(1:sim.net.n,1:sim.net.n,'k-');
set(gca,'xtick',linspace(1,sim.n,sim.n)-.5,'xticklabel',[],...
    'xgrid','on','xcolor','w',...
    'ytick',linspace(1,sim.n,sim.n)-.5,'ytickLabel',[],...
    'ygrid','on','ycolor','w',...
    'gridLineStyle','-','linewidth',2,'gridalpha',1)
xlabel('nodes'); ylabel('nodes');
figformat;
tl = title('Functional connections');tl.FontWeight = 'normal';

subplot(2,2,3);
mY       = squeeze(mean(sim.Y));
limY     = max(abs([min(mY(:)) max(mY(:))]));
limY     = limY+limY*.2;
plot(sim.time,mY,'linewidth',1.5);
ylim([-limY limY])
states   = unique(sim.summary.time);
states   = num2cell([states(1:end-1) states(2:end)],2);
col      = cbrewer('qual','Set1',max(numel(states),3));
col      = mat2cell(col(1:numel(states),:),ones(1,numel(states)),3);
hold on
ShadePlotForEmpahsis(states,col,0.1);
xlabel('time(s)/states'); ylabel('activity(a.u.)')
tl = title('Surrogate time-series'); tl.FontWeight = 'normal';
figformat % figformat(1,0,0.1);

[psd,f]  = multi_pwelch(mean(sim.Y),sim.srate);
subplot(2,2,4)
plot(f,squeeze(mean(pow2db(psd(:,:,:)),1)),'linewidth',1.5);hold off
figformat; %figformat(0,0,0.1);
xlabel('F(Hz)');ylabel('psd(db)')
tl = title({'Surrogate psd','(trial averaged)'});tl.FontWeight = 'normal';

