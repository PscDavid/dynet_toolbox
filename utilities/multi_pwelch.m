function [pxx,f] = multi_pwelch(matr,srate)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Power spectral density estimate of the input signal using Welch's 
% overlapped segment averaging estimator (input signal divided 
% into the longest possible segments to obtain as close to but not exceed 8
% segments with 50% overlap)
%
% Last update: 22.08.2019
%--------------------------------------------------------------------------
% INPUTs:
% - matr should be: [trialsXchannlesXtime]
% - srate:          sampling frequency
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fs = srate;
for tr = 1:size(matr,1)
    for ch = 1:size(matr,2)
        x       = squeeze(matr(tr,ch,:));
        [pxx(tr,ch,:),f] = pwelch(x',[],[],[],fs);
    end
end
