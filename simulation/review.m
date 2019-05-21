function review(sim)

figure('units','normalized','position',[0.25 0.1 .5 .8]); 
subplot(2,2,1)
imagesc(sim.SC);hold on
% plot(1:sim.net.n,1:sim.net.n,'k-');
set(gca,'xtick',linspace(1,sim.net.n,sim.net.n)-.5,'xticklabel',[],...
    'xgrid','on','xcolor','w',...
    'ytick',linspace(1,sim.net.n,sim.net.n)-.5,'ytickLabel',[],...
    'ygrid','on','ycolor','w',...
    'gridLineStyle','-','linewidth',2,'gridalpha',1)
xlabel('nodes');ylabel('nodes');
figformat;
tl = title('Structural links');tl.FontWeight = 'normal';
subplot(2,2,2)
imagesc(sim.DC);hold on
% plot(1:sim.net.n,1:sim.net.n,'k-');
set(gca,'xtick',linspace(1,sim.net.n,sim.net.n)-.5,'xticklabel',[],...
    'xgrid','on','xcolor','w',...
    'ytick',linspace(1,sim.net.n,sim.net.n)-.5,'ytickLabel',[],...
    'ygrid','on','ycolor','w',...
    'gridLineStyle','-','linewidth',2,'gridalpha',1)
xlabel('nodes');ylabel('nodes');
figformat;
tl = title('Functional channels');tl.FontWeight = 'normal';

subplot(2,2,3);
mY       = squeeze(mean(sim.Y));
limY     = max(abs([min(mY(:)) max(mY(:))]));
limY     = limY+limY*.2;
plot(sim.net.time,mY,'linewidth',1.5);
ylim([-limY limY])
states   = unique(sim.summary.time);
states   = num2cell([states(1:end-1) states(2:end)],2);
col      = cbrewer('qual','Set1',max(numel(states),3));
col      = mat2cell(col(1:numel(states),:),ones(1,numel(states)),3);
hold on
ShadePlotForEmpahsis(states,col,0.1);
xlabel('time(s)/states');ylabel('activity(a.u.)')
tl = title('Surrogate time-series');tl.FontWeight = 'normal';
figformat(1,0,0.1);

[psd,f]  = multi_pwelch(mean(sim.Y),sim.net.srate);
subplot(2,2,4)
plot(f,squeeze(mean(pow2db(psd(:,:,:)),1)),'linewidth',1.5);hold off
figformat(0,0,0.1);
xlabel('F(Hz)');ylabel('psd(db)')
tl = title({'Surrogate psd','(trial averaged)'});tl.FontWeight = 'normal';

