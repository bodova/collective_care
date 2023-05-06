function[] = plotrates(ratesA,typeModel,transitions,times,col,lw)
global idty 
global NP NR ND NL yP yR yD yL

onlyrates = 1;
switch typeModel
    case 9; xtitle = 'Transition'; ytitle = 'Rate';
    case 1; xtitle = 'P activity'; ytitle = 'Number'; ytitle2 = 'Rate';
    case 2; xtitle = 'R activity'; ytitle = 'Number'; ytitle2 = 'Rate';
    case 6; xtitle = 'load'; ytitle = 'Number'; ytitle2 = 'Rate';
    case 7; xtitle = 'log2(load)'; ytitle = 'P activity';
    case 8; xtitle = 'log2(load)'; ytitle = 'R activity';
end

titlN = {'Transitions $X \rightarrow S$','Transitions $X \rightarrow A$','Transitions $S \rightarrow X$','Transitions $S \rightarrow S$','Transitions $S \rightarrow A$','Transitions $A \rightarrow X$','Transitions $A \rightarrow S$','Transitions $A \rightarrow A$'};
titl = {'$X \rightarrow S$','$X \rightarrow A$','$S \rightarrow X$','$S \rightarrow S$','$S \rightarrow A$','$A \rightarrow X$','$A \rightarrow S$','$A \rightarrow A$'};
titl3 = {'$X$','$S$','$A$'};
%k = 3-idty; figure(k); %clf; % # transitions
figure(1); clf; hold on;
set(gcf,'color','w');

if typeModel == 9 %% constant
    if onlyrates == 1
        plot(exp(ratesA),'o','linewidth',1,'markersize',lw,'color',col,'markerfacecolor',col); hold on
        xticky = {'XS','XA','SX','SA','AX','AS','SS','AA'};
        xticks([1:8]);
        xticklabels(xticky);
        set(gca,'fontsize',14)
        xlabel(xtitle); ylabel(ytitle); xlim([1,8])
        x0=100; y0=100; width = 300; height = 250; set(gcf,'position',[x0,y0,width,height])
    else
        subplot(1,3,1); plot(transitions,'ok','linewidth',1,'markersize',lw,'markerfacecolor','k'); hold on
        set(gca,'fontsize',14); xlabel(xtitle); ylabel('Transitions (#)'); xlim([1,8])
        xticky = {'XS','XA','SX','SA','AX','AS','SS','AA'}; xticks([1:8]); xticklabels(xticky);

        subplot(1,3,2); plot(times/15,'ok','linewidth',1,'markersize',lw,'markerfacecolor','k'); hold on
        set(gca,'fontsize',14); xlabel('State'); ylabel('Time (s)'); xlim([1,3])
        xticky = {'X','S','A'}; xticks([1:3]); xticklabels(xticky);

        subplot(1,3,3); plot(exp(ratesA),'o','linewidth',1,'markersize',lw,'color',col,'markerfacecolor',col); hold on
        set(gca,'fontsize',14); xlabel(xtitle); ylabel(ytitle);
        xticky = {'XS','XA','SX','SA','AX','AS','SS','AA'}; xticks([1:8]); xticklabels(xticky); xlim([1,8])
        x0=100; y0=100; width = 1000; height = 250; set(gcf,'position',[x0,y0,width,height])
    end
end
    
if (typeModel == 1) || (typeModel == 2) || (typeModel == 4) || (typeModel == 6)
    if (typeModel == 1)  N = NP; y = yP; xlab = {'<0.01','medium','>0.2'};
    elseif (typeModel == 2)  N = NR; y = yR; xlab = {'<0.01','medium','>0.2'};
    elseif typeModel == 6 N = NL; y = yL; xlab = {'<5000','medium','>30000'};
    end
    
    id = [1,2,3,7,4,5,6,8];
    if idty == 2 cl = [0.5,0.5,0.5]; cl2 = [1,0,0];
    else cl = [0,0,0]; cl2 = [0.7,0.3,0.3];
    end
    
    for i=1:8
        %% rates
        D = exp(ratesA); Dmax = max(D(:));
        figure(1); subplot(3,3,i+1); hold on
        plot(1:3,D(id(i),:),'o-','color',col,'linewidth',1,'markersize',lw,'markerfacecolor',col);
        xticks(1:3); xticklabels(xlab)
        ylabel(ytitle2); xlabel(xtitle);
        ax = gca; ax.XTickLabelRotation = 0; ax.YTickLabelRotation = 0;
        title(titl(i),'Interpreter','latex'); box on;
        set(gca,'fontsize',14)
        x0=100; y0=100; width = 700; height = 500; set(gcf,'position',[x0,y0,width,height])

        %% transitions
        D = transitions; Dmax = max(D(:));
        figure(2); subplot(3,3,i+1); hold on
        plot(1:3,D(id(i),:),'o-','color','k','linewidth',1,'markersize',lw,'markerfacecolor','k');
        xticks(1:3); xticklabels(xlab)
        ylabel('Transitions (#)'); xlabel(xtitle);
        ax = gca; ax.XTickLabelRotation = 0; ax.YTickLabelRotation = 0;
        title(titl(i),'Interpreter','latex'); box on;
        set(gca,'fontsize',14)
        x0=100; y0=100; width = 700; height = 500; set(gcf,'position',[x0,y0,width,height])
    end
    for i=1:3
        %% times
        D = times; Dmax = max(D(:));
        figure(3); subplot(1,3,i); hold on
        plot(1:3,D(id(i),:),'o-','color','k','linewidth',1,'markersize',lw,'markerfacecolor','k');
        xticks(1:3); xticklabels(xlab)
        ylabel('Times (s)'); xlabel(xtitle);
        ax = gca; ax.XTickLabelRotation = 0; ax.YTickLabelRotation = 0;
        title(titl3(i),'Interpreter','latex'); box on;
        set(gca,'fontsize',14)
        x0=100; y0=100; width = 700; height = 200; set(gcf,'position',[x0,y0,width,height])
    end
    
elseif (typeModel == 7)|| (typeModel == 8) 
    if (typeModel == 7)  
        N1 = NP; N2 = NL; y1 = yP; y2 = yL; 
        xlab = {'<5000','medium','>30000'};
        ylab = {'<0.01','medium','>0.2'};
    elseif (typeModel == 8) 
        N1 = NR; N2 = NL; y1 = yR; y2 = yL;
        xlab = {'<5000','medium','>30000'};
        ylab = {'<0.01','medium','>0.2'};
    end
    
    if (typeModel == 1)  N = NP; y = yP; xlab = {'<0.01','medium','>0.2'};
    elseif (typeModel == 2)  N = NR; y = yR; xlab = {'<0.01','medium','>0.2'};
    elseif typeModel == 6 N = NL; y = yL; xlab = {'<5000','medium','>30000'};
    end

    id = [1,2,3,7,4,5,6,8];
    for i=1:8
        D = exp(ratesA);
        subplot(3,3,i+1);
        D = reshape((D(id(i),:,:)),[N1,N2]); Da = [D,zeros(N1,1)]; Db = [Da',zeros(N2+1,1)]';
        y1a = [2*y1(1)-y1(2),y1(1:end-1)]; y1b = (y1+y1a)/2; y1b = [y1b,y1(end)+(y1(end)-y1(end-1))/2];
        y2a = [2*y2(1)-y2(2),y2(1:end-1)]; y2b = (y2+y2a)/2; y2b = [y2b,y2(end)+(y2(end)-y2(end-1))/2];
        M = max(D(:));  pcolor(Db); 
        ColorDef = Col; colormap(ColorDef); colorbar;
        xlabel(xtitle); ylabel(ytitle);
        ylim([1,4])
        xlim([1,N2+1])
        title(titl(i),'Interpreter','latex');
        xticks(1.5:1:3.5); xticklabels(xlab)
        yticks(1.5:1:3.5); yticklabels(ylab)
    end
    x0=100; y0=100; width = 1000; height = 700; set(gcf,'position',[x0,y0,width,height])
end


