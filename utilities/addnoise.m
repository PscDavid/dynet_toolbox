function [Y,noise] = addnoise(Y,SNRdb,type)
%==========================================================================
% Add white or 1/f noise
%--------------------------------------------------------------------------
% [tr,m,~]    = size(Y);
Ps          = std(Y(:));
Pw          = Ps/10^((SNRdb/20));
switch type
    case 'w'
        noise        = Pw.*randn(size(Y));
        Y            = Y + noise;
    case '1/f'
        noise        = Pw.*randnd(-1,size(Y));
        Y            = Y + noise;
end            

