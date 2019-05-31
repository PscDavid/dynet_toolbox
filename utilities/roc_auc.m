function [auc,sens,spec] = roc_auc(X,Y,ncrit,figout)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Returns the Area-Under-the-ROC-Curve 
% and the sensitivity (true  positive rate)
% and specificity (1 - false positive rate)
% for a range of criteria based on quantile thresholding
%                                       D. Pascucci, University of Fribourg
% Last update: 19.05.2019
%--------------------------------------------------------------------------
% INPUT:
% -X:      ground truth
% -Y:      estimated connectivity matrix (e.g., PDC)
% -ncrit:  number of criteria used (spacing from 0-1)
% -figout: plot results
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

default('figout',0);
default('ncrit',20);

% check and remove diagonal elements
if ~isnan(Y(1)) && numel(size(Y))>2
    dg     = size(Y,1);
    for ij = 1:dg
        X(ij,ij,:,:) = NaN;
        Y(ij,ij,:,:) = NaN;
    end
end

rmv        = isnan(X);
X          = X(~rmv);
Y          = Y(~rmv);

qrange     = linspace(0.001,1,ncrit);
sens       = zeros(numel(qrange),1); spec    = sens;
for k = 1:numel(qrange)
    thre     = quantile(Y,qrange(k));
    sens(k,:)= mean((Y(X>0)>=thre));
    spec(k,:)= mean((Y(X==0)<thre));
end
auc        = trapz(1-sens, spec);

if figout
    plot(1-spec,sens,'-','linewidth',2);axis([0 1 0 1]);hold on
    d          = linspace(0,1,100);
    hold on;plot(d,d,'k--');
    xlabel('1-Specificity');
    ylabel('Sensitivity');
    set(gca,'xtick',[0 .5 1],'ytick',[0 .5 1]);
    tl  = title('ROC');tl.FontWeight = 'normal';
    figformat(0,0,.1);axis square
end

