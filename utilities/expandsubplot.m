function expandsubplot
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Expand (maximise) subplot figure temporarily ? then collapse it back 
% 
% Last update: 22.08.2019
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
set(gca, 'ButtonDownFcn',@(ax,~) ...
    set(copyobj(ax,uipanel('Position',[0,0,1,1])),...
    'Units','normal', 'OuterPosition',[0,0,1,1], ...
    'ButtonDownFcn',@(co,~)delete(get(co,'parent'))));
end