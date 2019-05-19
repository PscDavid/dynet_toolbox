
%%


%%
SM         = sim.SC;
FM         = sim.DC;
DM         = [randsample(1:10:200,sim.net.n)', ...
              randsample(1:10:200,sim.net.n)'];
[i,j]      = find(SM.*(eye(sim.net.n)==0));
plot(DM([i j],1),DM([i j],2),'k-','linewidth',2);hold on;
plot(DM(:,1),DM(:,2),'ok','linewidth',2,'markersize',40,...
                                               'markerfacecolor',[1 1 1]);
tx = regexprep(cellstr([repmat('n',[size(SM,1),1]) ...
                 num2str([1:size(SM,1)]')]),'\s','')';                                          
text(DM(:,1),DM(:,2),tx);

%%
SM         = sim.SC;
FM         = sim.DC;

clear i
theta = [0 : 2*pi/(sim.net.n) : 2*pi]
rad = 150
Pos = rad * exp(i*theta);
DM  = [real(Pos)' imag(Pos)'];DM(end,:)=[];
% plot the vertices on a circle
[jj,ii]      = find(SM.*(eye(sim.net.n)==0));
slink        = unique(sortrows([jj ii]),'rows');
for k = 1:numel(slink(:,1))
    plot(DM(slink(k,:),1),DM(slink(k,:),2),'k-','linewidth',2);hold on;
end
plot(DM(:,1),DM(:,2),'ok','linewidth',2,'markersize',25,...
                                         'markerfacecolor',[1 1 1]);
tx = regexprep(cellstr([repmat('n',[size(SM,1),1]) ...
                 num2str([1:size(SM,1)]')]),'\s','')';                                          
text(DM(:,1)-5,DM(:,2),tx);
axlim = [-rad rad]'+[-.2*rad .2*rad; -.2*rad .2*rad ]';
axis(axlim(:));
figformat;axis off

%%
clear i
theta = [0 : 2*pi/(sim.net.n) : 2*pi]
rad = 150
Pos = rad * exp(i*theta);
DM  = [real(Pos)' imag(Pos)'];DM(end,:)=[];
% plot the vertices on a circle
[jj,ii]      = find(SM.*(eye(sim.net.n)==0));
slink        = unique(sortrows([jj ii]),'rows');
for k = 1:numel(slink(:,1))
    X = DM(slink(k,:),1);
    Y = DM(slink(k,:),2);
    Xi = mean(X);
    Yi = mean(Y) + .2*mean(Y);  
    Xa = [X(1) Xi X(2)];
    Ya = [Y(1) Yi Y(2)];

    t  = 1:numel(Xa);
    ts = linspace(min(t),max(t),numel(Xa)*10); % has to be a fine grid
    xx = spline(t,Xa,ts);
    yy = spline(t,Ya,ts);

    plot(xx,yy,'r-','linewidth',2); hold on; % curve
    plot(xx(end)-10,yy(end),'r.','markersize',50);hold on
%     plot(X,Y,'or')        % end points
%     plot(Xi,Yi,'xr')      % intermediate point
hold on
end
plot(DM(:,1),DM(:,2),'ok','linewidth',2,'markersize',25,...
                                         'markerfacecolor',[1 1 1]);
tx = regexprep(cellstr([repmat('n',[size(SM,1),1]) ...
                 num2str([1:size(SM,1)]')]),'\s','')';                                          
text(DM(:,1)-5,DM(:,2),tx);
axlim = [-rad rad]'+[-.2*rad .2*rad; -.2*rad .2*rad ]';
axis(axlim(:));
figformat;axis off


% figure
% [i,j]      = find(FM.*(eye(sim.net.n)==0));
% plot(DM([i j],1),DM([i j],2),'k-','linewidth',2);hold on;
% plot(DM(:,1),DM(:,2),'ok','linewidth',2,'markersize',20,...
%                                                'markerfacecolor',[1 1 1]);
% 
% for k = 1:numel(i)
%     plot(DM([i(k) j(k)],1),DM([i(k) j(k)],2),'k-');hold on;
% end

%%
free
% global srate p time 
freqs      = (5:80)';
                 
% profile on
sim        = dynet_sim(5,[],2,15,.5,[],2,400);
% profile viewer

gt_PDC     = dynet_ar2pdc(sim,srate,freqs,'sPDC',[],[],1);

dynet_connplot(gt_PDC,time,freqs,[],[],[],sim.DC,1)

sim.summary

%
% PDC             = DYN_ar_to_pdc(sim,srate,freqs,'sPDC',0);
% conn_matr_plot(PDC,[],time,freqs,[.01 .99],jet,0,sim.DC,0)

%
% conn_matr_plot(gt_PDC,[],time,freqs,[.01 .99],jet,0,[],0)

KF         = dynet_SSM_KF(sim.Y,p,0.01);
kf_PDC     = dynet_ar2pdc(KF,srate,freqs,'sPDC',[],[],1);
dynet_connplot(kf_PDC,time,freqs,[],[],[],sim.DC,1)

% conn_matr_plot(kf_PDC,[],time,freqs,[.01 .99],jet,0,[],0)

lambda     = dynet_ISS_GCV(sim.Y,p,[],[],1);
SA         = dynet_SSM_SALK(sim.Y,p,lambda);
sa_PDC     = dynet_ar2pdc(SA,srate,freqs,'sPDC',[],[],1);
dynet_connplot(sa_PDC,time,freqs,[],[],[],sim.DC,1)

% conn_matr_plot(sa_PDC,[],time,freqs,[.01 .99],jet,0,[],0)

%%
close all
clear SPEC SENS
gt_PDC     = dynet_ar2pdc(sim,srate,freqs,'sPDC',[],[],0);
X          = gt_PDC(~isnan(gt_PDC));
KF         = dynet_SSM_KF(sim.Y,p,0.03);
kf_PDC     = dynet_ar2pdc(KF,srate,freqs,'sPDC',[],[],0);
Y          = kf_PDC(~isnan(gt_PDC));

qrange     = linspace(0.001,1,20);
for k = 1:numel(qrange)
    thre   = quantile(Y,qrange(k));
    SENS(k,:)= mean((Y(X>0)>=thre));
    SPEC(k,:)= mean((Y(X==0)<thre));
end
plot(1-SPEC,SENS,'-o');axis([0 1 0 1]);hold on
auc_kf = trapz(1-SENS, SPEC)


lambda     = dynet_ISS_GCV(sim.Y,p,[],[],0);
SA         = dynet_SSM_SALK(sim.Y,p,lambda);
sa_PDC     = dynet_ar2pdc(SA,srate,freqs,'sPDC',[],[],1);
Y          = sa_PDC(~isnan(gt_PDC));


[auc,sens,spec] = roc_auc(X,Y,50,1)

qrange     = linspace(0.001,1,20);
for k = 1:numel(qrange)
    thre   = quantile(Y,qrange(k));
    SENS(k,:)= mean((Y(X>0)>=thre));
    SPEC(k,:)= mean((Y(X==0)<thre));
end
plot(1-SPEC,SENS,'-o');axis([0 1 0 1]);hold on
auc_sa = trapz(1-SENS, SPEC)

d          = linspace(0,1,100);
hold on;plot(d,d,'k--');

%%
clist = logspace(-4,0,20);
qrange     = linspace(0.001,1,20);
for ck = 1:numel(clist)
    KF         = dynet_SSM_KF(sim.Y,p,clist(ck));
    kf_PDC     = dynet_ar2pdc(KF,srate,freqs,'sPDC',[],[],0);
    Y          = kf_PDC(~isnan(gt_PDC));

    for k = 1:numel(qrange)
        thre   = quantile(Y,qrange(k));
        SENS(k,:)= mean((Y(X>0)>=thre));
        SPEC(k,:)= mean((Y(X==0)<thre));
    end
    
    plot(1-SPEC,SENS,'-o');axis([0 1 0 1]);hold on
    drawnow
    auc(ck,:) = trapz(1-SENS, SPEC);
end
d          = linspace(0,1,100);
hold on;plot(d,d,'k--');



%%
free
scalefun   = @(x,y) (x-min(x(:)))/range(x(:))...
                     *range(y(:))+min(y(:));
                 

clist      = logspace(-3,log10(1),10);

n          = 5;      p          = 10;
srate      = 250;    duration   = 1.5;
sparsity   = 0.1;    nstates    = 3;
trials     = 300;

sim       = dynet_sim(n,srate,duration,p,sparsity,[],nstates,trials);

max_n      = 100;
noise      = normrnd(mean(sim.Y(:)),std(sim.Y(:)),...
                     trials,max_n-n,numel(sim.net.time));
                 
 
dim        = [5 7 10 20 30 50 100];%fix(logspace(log10(n),log10(max_n),5));  
RMSE       = NaN(numel(clist),numel(dim));
SPEC         = RMSE;
for d = 1:numel(dim)
    Y      = cat(2,sim.Y,noise(:,1:(dim(d)-n),:));
    lambda = dynet_ISS_GCV(Y,p,[],[],1);
    SA     = dynet_SSM_SALK(Y,p,lambda);    
    estAR  = SA.AR(1:n,1:n,:,:);
    outAR  = SA.AR((1+n):end,(1+n):end,:,:);
    lam(:,d) = lambda;
    SPEC_sa(:,d) = rms(outAR(:))
    RMSE_sa(:,d) = rms(estAR(:)-sim.AR(:));
    for k = 1:numel(clist)
        KF        = dynet_SSM_KF(Y,p,clist(k));
        estAR     = KF.AR(1:n,1:n,:,:);
        outAR     = KF.AR((1+n):end,(1+n):end,:,:);
        SPEC(k,d)   = rms(outAR(:));
        RMSE(k,d) = rms(estAR(:)-sim.AR(:));
    end
    prcdone(d,numel(dim),'c over n',20)
end

subplot(1,2,1)
plot(clist,RMSE);set(gca,'xscale','log')
subplot(1,2,2)
plot(clist,SPEC);set(gca,'xscale','log')