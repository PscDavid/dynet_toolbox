function dynet_connplot(ConnMatrix,time,freq,labels,quantrange,cmap,SC,univ)
% Beta version, needs to be commented and improved

% vline
% subtightplot

default('univ',0);   default('SC',0);
default('cmap',jet); default('quantrange',[0.01 0.99]);
default('labels',regexprep(cellstr([repmat('n',[size(ConnMatrix,1),1]) ...
                 num2str([1:size(ConnMatrix,1)]')]),'\s','')');

dim  = size(ConnMatrix);
if numel(quantrange)>1 && quantrange(1)>=0
maxscale = quantile(ConnMatrix(:),quantrange(2));
minval = min(ConnMatrix(:));
if minval<0
    minscale = quantile(ConnMatrix(:),quantrange(1));
else
    minscale = 0;
end
elseif numel(quantrange)>1 && quantrange(1)<0
    minscale = quantrange(1);
    maxscale = quantrange(2);
else
    minscale = 0;
    maxscale = quantrange; % here equal to an imposed maximum
end

figure('units','normalized','position',[0.25 0.2 .5 .6]); 

for rr = 1:dim(1)*dim(2)
    [i1, i2] = ind2sub([dim(1) dim(2)],rr);
    expandsubplot
    subtightplot(dim(1),dim(2),rr,[0.02 0.02],0.15,0.15)
%     subplot(dim(1),dim(1),rr)
    
    if univ==0
    if i1~=i2
       IM = squeeze(ConnMatrix(i2,i1,:,:));
       imagesc(time,freq,IM);colormap(cmap)
       axis xy
       vline(0,'w')       
    end
    else
        IM = squeeze(ConnMatrix(i2,i1,:,:));
        imagesc(time,freq,IM);colormap(cmap)
        axis xy
        vline(0,'w')       
    end
    
    if rr == dim(1)*dim(2) && univ==0
        axis off
        hp4 = get(subtightplot(dim(1),dim(2),rr),'Position');
        colorbar('Position', [hp4(1)+hp4(3)-.05   hp4(2)+0.1 ...
                                           0.02   hp4(2)+hp4(3)*1.5])
    end
    
    if i1 == 1 && i2 ~= 1
        ylabel(labels(i2))
        set(get(gca,'YLabel'),'Rotation',45);
    end
    if i2 == 1
        xlabel(labels(i1))    
        set(gca,'XAxisLocation','top');
        set(get(gca,'XLabel'),'Rotation',45);
    end
    if i2 == i1 && univ==0
        axis off
    end
    if rr == dim(1)
         ylabel(labels(1))    
         set(gca,'YAxisLocation','right')
    end
    if rr == (dim(1)*(dim(2)-1))+1 && univ==0
         xlabel(labels(1))    
         set(gca,'XAxisLocation','bottom')
    end
    
    if rr == (dim(1)*(dim(2)-1)) && univ==0
        axis on
        xlabel('time (s)');grid on
        ylabel('f (Hz)')
    else
        set(gca,'YTick',[],'XTick',[])       
    end
    
    if univ && rr==1
        axis on
        set(gca,'XTick',[]) 
        yticks('auto');grid on
        ylabel('f (Hz)')
    end
    if univ && rr==(dim(1)*dim(2))
        axis on
        set(gca,'YTick',[]) 
        xticks('auto');grid on
        xlabel('time (s)')
        hp4 = get(gca,'Position');
        colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2) ...
                                            0.02  hp4(2)+hp4(3)*1.5])
    end
        
%     set(gca,'FontSize',12)
    figformat(0,0,.1,10);
    if nnz(SC)>0 && SC(i2,i1)==1 && i1~=i2
        ax = gca;
        ax.XColor = [0.8500 0.3250 0.0980];%[0.9290 0.6940 0.1250]; 
        ax.YColor = [0.8500 0.3250 0.0980];%[0.9290 0.6940 0.1250];
        ax.LineWidth = 2;
        xc = get(ax,'xlabel');        xc.Color  = [0 0 0];
        yc = get(ax,'ylabel');        yc.Color  = [0 0 0];
    else
        ax = gca;
        ax.XColor = [0 0 0];%[0.8500 0.3250 0.0980];%[0.9290 0.6940 0.1250]; 
        ax.YColor = [0 0 0];%[0.8500 0.3250 0.0980];%[0.9290 0.6940 0.1250];
        ax.LineWidth = 2;
    end        
        
%     colormap(cmap)
    caxis([minscale maxscale])
end
% axis off

