function [auc,sens,spec] = roc_auc(X,Y,ncrit,figout)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Returns the Area-Under-the-ROC-Curve 
% and the sensitivity (true  positive rate; X,Y>0)
% and specificity (1 - false positive rate; X,Y==0)
% for a range of criteria based on quantile thresholding
%                                       
% Last update: 22.08.2019
%--------------------------------------------------------------------------
% INPUTs:
% -X:      ground truth
% -Y:      estimated connectivity matrix (e.g., PDC)
% -ncrit:  number of criteria used (spacing from 0-1)
% -figout: defaul 0 othersiwe 1 to plot results
%--------------------------------------------------------------------------
% INVOKED FUNCTIONs
% default.m; figformat.m
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
sens       = zeros(numel(qrange),1); 
spec       = sens;
thre       = quantile(Y,qrange);

for k = 1:numel(qrange)
    sens(k,:)= mean((Y(find(X>0))>=thre(k)));
    spec(k,:)= mean((Y(find(X==0))<thre(k)));
end

auc        = abs(trapz(1-spec, sens));

if figout
    plot(1-spec,sens,'-','linewidth',2)
    axis([0 1 0 1])
    hold on
    d      = linspace(0,1,100);
    plot(d,d,'k--');
    xlabel('1-Specificity');
    ylabel('Sensitivity');
    set(gca,'xtick',[0 .5 1],'ytick',[0 .5 1]);
    tl  = title('ROC');
    tl.FontWeight = 'normal';
    figformat(0,0,.1);
    axis square
end

