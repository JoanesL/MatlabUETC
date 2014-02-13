function PlotLogNegative(x,y,colour,lineWidth,minLevel)

if exist('lineWidth','var')~=1; lineWidth=1; end

if nargin<3
    colour='b'
end
if exist('minLevel','var')~=1; minLevel=1e-5*max(abs(y)); end

if min(y)<=0
    yPlot=y(1);
    xPlot=x(1);
    for i=2:max(size(y))
        if sign(y(i-1))~=sign(y(i))
            yPlot=[yPlot 0 y(i)];
            xPlot=[xPlot ( x(i-1)*abs(y(i-1))+x(i)*abs(y(i)) ) / ( abs(y(i-1))+abs(y(i)) ) x(i)];
        else
            yPlot=[yPlot y(i)];
            xPlot=[xPlot x(i)];
        end
    end
    
    yPos=yPlot; yPos(yPos<=0)=minLevel;
    plot(xPlot,yPos,['.-' colour],'LineWidth',lineWidth); hold on
    yNeg=-yPlot; yNeg(yNeg<=0)=minLevel;
    plot(xPlot,yNeg,['.--' colour],'LineWidth',lineWidth);
else
    plot(x,y,['.-' colour],'LineWidth',lineWidth)
end

disp(['maximoa ' int2str(max(abs(y))) ])

set(gca,'YScale','log')
set(gca,'YLim',[minLevel*1.001 max(abs(y))])
