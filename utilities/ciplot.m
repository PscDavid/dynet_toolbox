function hp = ciplot(x,y,color,dots)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Figure displaying y with shadow confidence intervals and distribution
%
% Last update: 22.08.2019
%--------------------------------------------------------------------------
% INPUTs
% - x:     x-axis vector
% - y:     y-axis vector
% - color: color specification (see ColorSpec Mathworks)
% - dots:  0, default otherwisw 1 for log-space x
%--------------------------------------------------------------------------
% INVOKED FUNCTIONs
% default.m
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
default('dots',0);

my              = mean(y);
bci             = bootci(10000,{@mean,y});
plot(x,bci','color',color)
hold on
% plot(x,my,'-o','color',color);hold on
hp              = patch([x'; x(end:-1:1)'; x(1)], [bci(1,:)'; ...
    bci(2,end:-1:1)'; bci(1)], color);
hp.FaceAlpha    = 0.5;
% hold on;
hl              = line(x,my);
hl.Color        = color;
hl.LineWidth    = 1.5;
if dots % beta option, for log spaced x
    jitx            = x+normrnd(x,x.*.1,1,numel(x)).*...
        normrnd(0,.01,size(y,1),numel(x));
    szy             = 1./abs(y-my);
    szy             = (50*(szy./max(szy)));
    for k = 1:numel(x)
        scatter(jitx(:,k),y(:,k),szy(:,k),color)
        hold on
    end
end
