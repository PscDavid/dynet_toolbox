free
% network properties
n        = 10;
p        = 7;
snr      = 1; % write NaN, if only signal
sim      = dynet_sim(n,[],[],p,[],[],[],snr);
popt     = size(sim.AR,3);
Y        = sim.Y;
% search grid for p
listap   = 1:1:20;
pVar     = 0.99; % variance criterion in SALKt
c        = 0.01;% adaptation constant in KF
for k = 1:numel(listap)
    p      = listap(k);
    SA     = dynet_SSM_SALKff(Y,p,pVar);
    KF     = dynet_SSM_KF(Y,p,c);
    
    tmpY   = Y;
    tmpY(:,:,1:p) = 0; % remove p-lags from original signals
    RMS_SA(k,:) = rms(tmpY(:)-SA.PY(:));
    RMS_KF(k,:) = rms(tmpY(:)-KF.PY(:));
    disp(['p ' num2str(p) ' of max p= ' num2str(listap(end))]);
end

plot(listap,RMS_KF,'ko-');hold on
plot(listap,RMS_SA,'ro-');
xlabel('p');ylabel('rms');
legend('KF','SALKt');
vline(popt,'b--')
figformat