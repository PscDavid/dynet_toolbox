function lambda_gcv = dynet_ISS_GCV(Y,p,nISS,lvec,figout)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Independent Segment Sampling for lambda selection through the
% Generalized cross-validation method [1].
%                                       D. Pascucci, University of Fribourg
% Last update: 17.05.2019
% INPUT:
% Y    = data                (trials x nodes x time);
% nISS = 100;                number of independent segments
% p    = 5;                  model order
% lvec = logspace(-4,0,50);  lambda grid (percentage of infinity norm)
%--------------------------------------------------------------------------
% [1] Golub, G. H., Heath, M., & Wahba, G. (1979). 
% Generalized Cross-Validation as a Method for Choosing a Good Ridge
% Parameter. Technometrics, 21(2), 215–223.

%--------------------------------------------------------------------------
% default
default('nISS',50);
default('lvec',logspace(-4,0,20));
default('figout',0);

%--------------------------------------------------------------------------
% define segments of data
[trl,dim,tm] = size(Y);
v            = ((p+1):tm)';
xr           = randsample(v,nISS)+(-1:-1:-p); %predictor segments(past Y)
yr           = xr(:,1)+1;%index for actual data
I            = eye(trl);
GCV          = NaN(numel(yr),numel(lvec));

%--------------------------------------------------------------------------
% loop over segments
for j = 1:numel(yr)
    % get segments
    segY   = Y(:,:,yr(j,:));
    segX   = reshape(Y(:,:,xr(j,:)),[trl dim*p]);
    for k = 1:numel(lvec)
        pc_lambda  = lvec(k);
        lambda_max = 0.8*norm(segX'*segY,'inf'); % if lambda_max = norm(segX'*segY,'inf')-> solution=0
        lambda_k   = pc_lambda*lambda_max;
        % economy-size decomposition of trl-by-dim H, 
        [U,S,V]    = svd(segX,'econ');
        % only the first dim columns of U are computed, and S is dim-by-dim
        d          = diag(S);
        D          = diag(d./(d.^2+lambda_k));
        Hinv       = V*D*U';
        % hat/influence matrix
        % Sherman-Morrison-Woodbury formula 
        A          = segX*Hinv; %X*pinv(X'*X + lambda_k*I)*X';
        den        = (1/(trl))*trace(I-A).^2;
        RSS        = mean(sum(((I-A)*segY).^2)); % RSS
        GCV(j,k)   = RSS / den;
    end
end
mGCV        = mean(GCV);
min_lambda  = lvec(dsearchn(mGCV',min(mGCV)));
from_minGCV = mGCV(dsearchn(mGCV',min(mGCV)):end);
idx_opt_lam = dsearchn(mGCV',min(mGCV))+...
              (dsearchn(from_minGCV',min(mGCV)+0.5*std(mGCV))-1);
lambda_gcv  = lvec(idx_opt_lam);

if figout
    figure
    plot(lvec,mean(GCV),'k-','linewidth',1.5);
    set(gca,'xscale','log');
    xlabel('\lambda');ylabel('GCV')
    hold on;
    pl1=plot(min_lambda,mGCV(dsearchn(mGCV',min(mGCV))),'.k','markersize',40);
    pl2=plot(lambda_gcv,mGCV(idx_opt_lam),'.r','markersize',40);
    legend([pl1 pl2],{'min','0.5 std'})
    figformat
end