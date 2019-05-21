function [pxx,f] = multi_pwelch(matr,srate)
% matr should be: trialsXchXtime
% before was chXtimeXtrials
fs = srate;
for tr = 1:size(matr,1)
    for ch = 1:size(matr,2)
        x       = squeeze(matr(tr,ch,:));        
        [pxx(tr,ch,:),f] = pwelch(x',[],[],[],fs);   
        % By default, x is divided into the longest possible segments to 
        % obtain as close to but not exceed 8 segments with 50% overlap. 
    end
end
% plot(f,pow2db(pxx))