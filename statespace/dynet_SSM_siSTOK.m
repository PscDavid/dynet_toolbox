function STOK = dynet_SSM_siSTOK(Y,p,C,ff,scale_f)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Self-Tuning Optimized Kalman filter with connectivity priors,
% Generalized Tikhonov Regularization
%                                                         D. Pascucci, EPFL
% Last update: 20.02.2021
%--------------------------------------------------------------------------
% INPUT:
% -Y:      matrix
%           [Trials X Channels X Time-samples]
%           Sample data
% -p:      scalar positive
%           p-order
% -C:      weighted/binary SC priors
% -ff:     TSVD filtering factor [2]
%--------------------------------------------------------------------------
% OUTPUT:   STOK, structure with fields:
% -AR:     matrix
%           [Channels X Channels X p X Time-samples]
%           MVAR model coefficients
% -R:      Measurement Noise Covariance Matrix (innovation based)
%           [Channels X Channels X Time]
% -PY:     matrix
%           [Trials X Channels X Time-samples]
%           One-step Prediction on Y
% -c:      self-tuning c for each k
%==========================================================================
% References:
% [1] Nilsson, M. (2006). Kalman filtering with unknown noise covariances.
%     In Reglermöte 2006.
% [2] Hansen, P. C. (1987). The truncated svd as a method for 
%     regularization. BIT Numerical Mathematics, 27(4), 534-553.
% [3] Yan, X., & Su, X. (2009). Linear regression analysis: theory 
%     and computing. World Scientific. (section 9.1.3)
% [4] Kaipio, J., & Somersalo, E. (2006). Statistical and computational 
%     inverse problems (Vol. 160). Springer Science & Business Media.
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -check input
default('ff',0.99); 
default('scale_f',[10^-4 0.1])
assert(numel(size(Y))==3,'Check the dimensions of Y');
[trl,dim,tm] = size(Y);
% -Preallocating main variable
AR           = zeros(dim,dim*p,tm);          % AR matrix
R            = zeros(dim,dim,tm);            % Additive observation noise
PY           = zeros(size(Y));               % One-step predictioR
% -Initializing variables
xm           = zeros(dim*p,dim)+1e-4;        % Prior state estimates
trEk         = zeros(tm,1);                  % Innovation monitoring
allc         = zeros(tm,1);                  % Self-tuning c, for all k
Cp           = repmat(C,[1 1 p]);            % Prior matrix
Cp           = rescale(Cp(:,:)',scale_f(1),scale_f(2));  % 10^-4 - 0.1
Q            = cell(1,dim);
for j = 1:dim
    Q{j}     = pinv(diag(Cp(:,j)));
end
%% ========================================================================
% -Loop from p(lag)+1 over time
for k = (p+1):tm
    %----------------------------------------------------------------------
    % Transition matrix H at time k
    data_back  = squeeze(Y(:,:,k-1:-1:k-p));
    H          = reshape(data_back,trl,dim*p);    
    % Measurement and Innovation
    Z          = Y(:,:,k);
    PY(:,:,k)  = H*xm;
    vk         = Z - PY(:,:,k); % Innovation at time k
    % Measurements innovation monitoring
    tmp        = (vk'*vk);
    R(:,:,k)   = tmp./max(trl-1,1);
    trEk(k,:)  = trace(tmp);
    % Damped SVD denoising 
    [H,~]      = tsvd_reg(H,ff);
    %% --------------------------------------------------------------------
    % Structural priors
    HH         = H'*H;
    HZ         = mat2cell(H'*Z,dim*p,ones(1,dim));
    betai_j    = cellfun(@(x,y)inverse(HH+x)*y,Q,HZ,'UniformOutput',false);
    betas      = cell2mat(betai_j);
    %% --------------------------------------------------------------------
    % Self-tuning adaptation constant
    if k>(p+1)*2 % 0.05<=c<=0.95
        ntrEk  = trEk(k-1:-1:(k-p*2));
        e_k    = mean(ntrEk(1:p));
        e_p    = mean(ntrEk(p+1:end));
        x      = abs(e_k-e_p)./e_p;
        c      = min(0.05 + x ,0.95); % or  x./max percentage allowed
    else
        c      = 0.05;
    end
    allc(k,:)  = c;
    %----------------------------------------------------------------------
    % Kalman update
    xp         = (xm + c*betas) ./ (1+c);  % [1]
    xm         = xp;
    AR(:,:,k)  = xm';
       
end
%==========================================================================
% -Saving output variables
AR             = reshape(AR,[dim dim p tm]);
STOK.AR        = AR;
STOK.R         = R;
STOK.PY        = PY;
STOK.c         = allc;
end

function [H,iH]  = tsvd_reg(H,ff)
    % Regularizing H by TSVD spectral smoothing [2]
    [U,S,V]    = svd(H,'econ');
    d          = diag(S);
    % fixed filtering factor threshold
    relv       = d.^2./sum(d.^2);
    filtfact   = find(double(cumsum(relv)<ff),1,'last');
    if isempty(filtfact)
        filtfact = 1;
    end
    lambda_k   = d(filtfact).^2;
    %diag(1./d.*((d.^2)./(d.^2+lambda_k)));
    D          = diag(d./(d.^2+lambda_k)); 
    iD         = diag(1./diag(D));%inverse(D);
    H          = U*iD*V'; % reconstruct H after SVD filtering
    if nargout==2
        iH     = V*D*U';
    end
end