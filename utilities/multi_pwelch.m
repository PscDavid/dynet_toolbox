function [pxx,f] = multi_pwelch(matr,srate)
% matr should be: trialsXchXtime
% before was chXtimeXtrials
fs = srate;
for tr = 1:size(matr,1)
    for ch = 1:size(matr,2)
        x       = squeeze(matr(tr,ch,:));        
        [pxx(tr,ch,:),f] = pwelch(x',[],[],[],fs);        
    end
end
% plot(f,pow2db(pxx))