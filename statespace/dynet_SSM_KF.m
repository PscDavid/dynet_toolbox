function KF = dynet_SSM_KF(Y,p,uc)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Kalman filter for state-space modeling of physiological time series
%                                       D. Pascucci, University of Fribourg
% Last update: 20.09.2018
%--------------------------------------------------------------------------
% INPUT:
% -Y:      matrix
%           [Trials X Channels X Time-samples]
%           Sample data
%           Classical for Trials = 1; General for Trials > 1 [1]
% -p:      scalar positive
%           p-order
% -c:      update Constants 0<=c<=1
%           [c C], C=c if numel(uc)==1
%--------------------------------------------------------------------------
% OUTPUT:   KF, structure with fields:
% -AR:     matrix
%           [Channels X Channels X p X Time-samples]
%           MVAR model coefficients
% -NS:     Measurement Noise Covariance Matrix (W in [1])
%           [Channels X Channels X Time]
% -PY:     matrix
%           [Trials X Channels X Time-samples]
%           One-step Predictions on Y
% -PYe:    matrix
%           [Trials X Channels X Time-samples]
%           One-step Residuals
% -c:      self-tuning c for each k
%==========================================================================
% References:
%
% [1] Milde, T., Leistritz, L., ..., & Witte, H. (2010). 
%     A new Kalman filter approach for the estimation of high-dimensional
%     time-variant multivariate AR models and its application in analysis
%     of laser-evoked brain potentials. Neuroimage, 50(3), 960-969.
% [2] Gelb, A. (Ed.). (1974). Applied optimal estimation. MIT press.
%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Notes: all equations are based on [1], but using the standard variable
% naming in Kalman filters [2]
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -check input
if numel(size(Y))<3;error('Check the dimensions of Y');end
if numel(uc)==1; C = uc; c = uc;
else; C = uc(1); c = uc(2); end

% -Preallocate main variable
[trl,dim,tm] = size(Y);
AR           = zeros(dim,dim*p,tm);          % AR matrix
R            = zeros(dim,dim,tm);            % Estimated noise measurement 
PY           = zeros(size(Y));               % One-step predictions

% -Initialize unknown design variables
Rk           = eye(dim)*1e-2;                % Noise measurement matrix
xm           = zeros(dim*p,dim)+1e-4;        % Prior state estimates
Pmi          = eye(dim*p)*1e-4;              % Estimated prediction error 
%==========================================================================
% -RECURSION
%   loop from p(lag)+1 through time
for k = (p+1):tm
    
    % Select previous data for the matrix H        
    data_back    = squeeze(Y(:,:,k-1:-1:k-p));
    H            = reshape(data_back,trl,dim*p);                           % see eq. (6)
    Z            = Y(:,:,k);      % current observation                    % see eq. (5)
    
    % Recursion on measurement noise covariance W
    Rn           = (Z-H*xm)'*(Z-H*xm)/max(trl-1,1);    
    Rk           = Rk+c*(Rn-Rk);  % adaptation constant on the past of R   % see eq. (13)
  
    % One-Step Prediction 
    PY(:,:,k)    = H*xm;
    
    % Innovation or residual covariance
    S            = H*Pmi*H' + sum(diag(Rk))*eye(trl);                      % see eq. (6b)
    
    % Optimal Kalman Gain 
    K            = Pmi*H'/S;                                               % see eq. (7)

    % A posteriori state estimate
    xp           = xm + K*(Z-H*xm);                                        % see eq. (8)
    xm           = xp;          % update xm                                % see eq. (9)
    
    % A posteriori estimate prediction error covariance
    Ppl         = Pmi-K*S*K';                                              % see eq.(10) 
    Pmi         = Ppl;          % update Pmi (a priori P)       
    Pmi(1:dim*p+1:end) = Pmi(1:dim*p+1:end) + C^2;                         % see eq.(11) 
   
    % Collect recursive estimates
    AR(:,:,k)   = xm';
    R(:,:,k)    = Rk;
    
end

% -Saving output variables
AR           = reshape(AR,[dim dim p tm]);
KF.AR        = AR;
KF.R         = R;
KF.PY        = PY;

