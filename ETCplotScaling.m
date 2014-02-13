%ETC plotting: produces figures like fig 8 of 2010 paper

function ETCplotScaling(All,id,run,tRef,tOffSet,xiscaling,tLimit,inPath,fitlimit,parent)

if nargin==0; 
  help ETCplot
  return
end

global gpath

if ~exist('inPath','var'); inPath=''; end 

if prod(size(inPath))>0; 
  path=inPath; 
else
  if prod(size(gpath))>0
    path=gpath;
  else
    disp(['Please set gpath global variable to default path'...
	  ' or specify path in fucntion call'])
    return
  end
end

if ~exist('tOffSet','var'); tOffSet=0; end
if ~exist('tLimit','var'); tLimit=[0 9999999]; end

%Prepare for plot
%if exist('parent','var')~=1; clf; else axes(parent); end

%Get number of runs
nRuns=size(run,2);

%Get tOffset from statsFile if necessary
if strcmp(tOffSet,'*')==1
  disp(['** Getting tOffSet from statsFile Lag. fit for ' ...
	'tRef -> 2*tRef **'])
  tOffSet = statsFile(-1,id,run,tRef*[1 2]);
end

if (xiscaling == 1)
    disp(['** Scaling with xiLag'])
    [xiLag tStat] = statGet('xiLag',id,run,path);
    if nRuns > 1
        xiLagAv = mean(xiLag,1);
    else
        xiLagAv = xiLag;
    end
end

%Get xiLag from statsFile 
%if strcmp(tOffSet,'xiscaling')==1
  
%  xiscaling = 1;
%end

%Duplicate tOffSet if single value (eg. 0) given for many runs
if size(tOffSet,2)==1 && nRuns>1
  tOffSet = ones(1,nRuns)*tOffSet;
end

%Plot all of them of top of each other

if All==1  
          
    for i=1:5
        %Load ETC
        if i==1
            Cname='scalar11';
        elseif i==2
            Cname='scalar12';
        elseif i==3
            Cname='scalar22';
        elseif i==4
            Cname='vector';
        elseif i==5
            Cname='tensor';
        end

    [k,t,C,sd]=ETCload(path,Cname,id,run,tRef,tOffSet,tLimit(1),tLimit(2));
    C=abs(C);

    if exist('number','var')==1
     C=C(number(1):number(2),:);
     t=t(number(1):number(2));
    end
    
    if (xiscaling == 1)
    %This must be done since usually tStat(end)<t(end), therefore it is
    %imposible to perform the interpolation to the lastest times.
        which=find(t<tStat(end));
        t=t(which);
        C=C(which,:);
        for i=1:size(C,1)
            xiScale = interp1(tStat,xiLagAv,t(i));
            C(i,:) = xiScale*C(i,:)/t(i);
            t(i) = xiScale;
        end
    end

    kt=k*t(1);
    zein = find( kt>=fitlimit(1) & kt<=fitlimit(2) );

    if strcmp(Cname,'vector')==1
    plotLogNegative(k*t(1),((k*t(1))).*C(1,:),'m',2)
    hold on;
    
    P = polyfit(log10(kt(zein)),log10(kt(zein).*C(1,zein)),1);
    Fit = P(1)*log10(kt) + P(2);

    plotLogNegative(kt,10.^(Fit),'-k',0.5);
    end
    if strcmp(Cname,'tensor')==1
    plotLogNegative(k*t(1),(k*t(1)).*C(1,:),'b',2)
    hold on;
    P = polyfit(log10(kt(zein)),log10(kt(zein).*C(1,zein)),1);
    Fit = P(1)*log10(kt) + P(2);

    plotLogNegative(kt,10.^(Fit),'-k',0.5);
    end
    if strcmp(Cname,'scalar12')==1
    plotLogNegative(k*t(1),(k*t(1)).*C(1,:),'r',2)
    hold on;
    P = polyfit(log10(kt(zein)),log10(kt(zein).*C(1,zein)),1);
    Fit = P(1)*log10(kt) + P(2);

    plotLogNegative(kt,10.^(Fit),'-k',0.5);
    end
    if strcmp(Cname,'scalar22')==1
    plotLogNegative(k*t(1),(k*t(1)).*C(1,:),'g',2)
    hold on;
    P = polyfit(log10(kt(zein)),log10(kt(zein).*C(1,zein)),1);
    Fit = P(1)*log10(kt) + P(2);

    plotLogNegative(kt,10.^(Fit),'-k',0.5);
    end
    if strcmp(Cname,'scalar11')==1
    plotLogNegative(k*t(1),(k*t(1)).*C(1,:),'k',2)
    hold on;
    P = polyfit(log10(kt(zein)),log10(kt(zein).*C(1,zein)),1);
    Fit = P(1)*log10(kt) + P(2);

    plotLogNegative(kt,10.^(Fit),'-k',0.5);
    end
%end
hold on
    x1=linspace(fitlimit(1),fitlimit(1),9000);
    x2=linspace(fitlimit(2),fitlimit(2),9000);
    y1=logspace(-5.5,4,9000);
    
    plot(x1,y1,'k'); hold on;
    plot(x2,y1,'k'); hold on;
        
%Highlight lines and plot uncertainties
%plotLogNegative(k*t(1),(k*t(1))*C(1,:),'b',2) 
%if nRuns>1 
%    plotErrorBars(k*t(1),(k*t(1))*C(1,:),sd(1,:),'b')
%end

set(gca,'XLim',[3 350])
set(gca,'YLim',[3 400])
%axis tight
set(gca,'XScale','log')
set(gca,'YScale','log')

if(xiscaling==1)
    xlabel('k\xi')
    ylabel('k\xi*C(k\xi,k\xi) or (k\xi)^3*C(k\xi,k\xi)')
else
    xlabel('kt')
    ylabel('kt*C(kt,kt) or (kt)^3*C(kt,kt)')
end

end

%EVOLUTION OF SLOPES
    figure()
    ax=multiPlot([2 3]);
    for i=1:5
        
        if i==1
        Cname='scalar11';
        elseif i==2
        Cname='scalar12';
        elseif i==3
        Cname='scalar22';
        elseif i==4
        Cname='vector';
        elseif i==5
        Cname='tensor';
        end
        
        [k,t,C,sd]=ETCload(path,Cname,id,run,tRef,tOffSet,tLimit(1),tLimit(2));
        C=abs(C);

        if exist('number','var')==1
         C=C(number(1):number(2),:);
         t=t(number(1):number(2));
        end
        
        if (xiscaling == 1)
        %This must be done since usually tStat(end)<t(end), therefore it is
        %imposible to perform the interpolation to the lastest times.
            which=find(t<900);
            t=t(which);
            C=C(which,:);
            for l=1:size(C,1)
                xiScale = interp1(tStat,xiLagAv,t(l));
                C(l,:) = xiScale*C(l,:)/t(l);
                time(l)=t(l);
                t(l) = xiScale;
            end
        end
        %HEMEN ATERA MALDA GUZTIAK
    
        for j=1:size(t,1)
            kt=k*t(j);
            zein = find( kt>=fitlimit(1) & kt<=fitlimit(2) );
            P = polyfit(log10(kt(zein)),log10(kt(zein).*C(j,zein)),1);
            Fit = P(1)*log10(kt) + P(2);
        
            %Store slopes
    
            x(i,j)=P(1);
        end
    end
        
    for j=1:5
    %Fit the slope
    s = fitoptions('Method','NonlinearLeastSquares');
    f = fittype('a + b/x^c','options',s);
    %[c2,gof2] = fit(t(:),(x(j,:))',f)

        if j==1
        axes(ax(1))        
        plot(time(:),x(j,:),'k'); hold on;
        %plot(c2,'-k')
        %hleg = legend(num2str(coeffvalues(c2)),'Location','SouthEast')
        % Make the text of the legend italic and color it brown
        %set(hleg,'FontAngle','italic','TextColor',[.3, .2, .1])
        %set(gca,'YLim',[-0.001 0.001])
        elseif j==2
        axes(ax(2))    
        plot(time(:),x(j,:),'r'); hold on;
        %plot(c2,'-k')
        %hleg = legend(num2str(coeffvalues(c2)),'Location','SouthEast')
        % Make the text of the legend italic and color it brown
        %set(hleg,'FontAngle','italic','TextColor',[.3, .2, .1])
        %set(gca,'YLim',[-0.001 0.001])
        elseif j==3
        axes(ax(3))
        plot(time(:),x(j,:),'g'); hold on; 
        %plot(c2,'-k')
        %hleg = legend(num2str(coeffvalues(c2)),'Location','SouthEast')
        % Make the text of the legend italic and color it brown
        %set(hleg,'FontAngle','italic','TextColor',[.3, .2, .1])
        %set(gca,'YLim',[-0.001 0.001])
        elseif j==4
        axes(ax(4))
        plot(time(:),x(j,:),'m'); hold on;
        %plot(c2,'-k')
        %hleg = legend(num2str(coeffvalues(c2)),'Location','SouthEast')
        % Make the text of the legend italic and color it brown
        %set(hleg,'FontAngle','italic','TextColor',[.3, .2, .1])
        %set(gca,'YLim',[-0.001 0.001])
        elseif j==5
        axes(ax(5))
        plot(time(:),x(j,:),'b'); hold on; 
        %plot(c2,'-k')
        %hleg = legend(num2str(coeffvalues(c2)),'Location','SouthEast')
        % Make the text of the legend italic and color it brown
        %set(hleg,'FontAngle','italic','TextColor',[.3, .2, .1])
        %set(gca,'YLim',[-0.02 -0.01])

        end
    end
    multiPlotZoom(ax);



%Plot Individually
else
    
    for i=1:5
    %Load ETC
if i==1
    Cname='scalar11';
elseif i==2
        Cname='scalar12';
elseif i==3
        Cname='scalar22';
elseif i==4
    Cname='vector';
elseif i==5
    Cname='tensor';
end


[k,t,C,sd]=ETCload(path,Cname,id,run,tRef,tOffSet,tLimit(1),tLimit(2));
C=abs(C);

    if exist('number','var')==1
     C=C(number(1):number(2),:);
    t=t(number(1):number(2));
    end

    
    %Select range to fit
    kt=k*t(1);
    zein = find( kt>=fitlimit(1) & kt<=fitlimit(2) );
   
    figure();

    if strcmp(Cname,'vector')==1
    plotLogNegative(k*t(1),((k*t(1))).*C(1,:),'k',2)
    hold on;
    
    P = polyfit(log10(kt(zein)),log10(kt(zein).*C(1,zein)),1);
    Fit = P(1)*log10(kt) + P(2);

    plotLogNegative(kt,10.^(Fit),'-k',0.5);
    end
    
    if strcmp(Cname,'tensor')==1
    plotLogNegative(k*t(1),(k*t(1)).*C(1,:),'b',2)
    hold on;
    P = polyfit(log10(kt(zein)),log10(kt(zein).*C(1,zein)),1);
    Fit = P(1)*log10(kt) + P(2);

    plotLogNegative(kt,10.^(Fit),'-k',0.5);
    end
    
    if strcmp(Cname,'scalar12')==1
    plotLogNegative(k*t(1),(k*t(1)).*C(1,:),'r',2)
    hold on;
    P = polyfit(log10(kt(zein)),log10(kt(zein).*C(1,zein)),1);
    Fit = P(1)*log10(kt) + P(2);

    plotLogNegative(kt,10.^(Fit),'-k',0.5);
    end
    
    if strcmp(Cname,'scalar22')==1
    plotLogNegative(k*t(1),(k*t(1)).*C(1,:),'g',2)
    hold on;
    P = polyfit(log10(kt(zein)),log10(kt(zein).*C(1,zein)),1);
    Fit = P(1)*log10(kt) + P(2);

    plotLogNegative(kt,10.^(Fit),'-k',0.5);
    end
    
    if strcmp(Cname,'scalar11')==1
    plotLogNegative(k*t(1),(k*t(1)).*C(1,:),'k',2)
    hold on;
    P = polyfit(log10(kt(zein)),log10(kt(zein).*C(1,zein)),1);
    Fit = P(1)*log10(kt) + P(2);

    plotLogNegative(kt,10.^(Fit),'-k',0.5);
    end
    
    hold on
    
    %Vertical lines
    x1=linspace(fitlimit(1),fitlimit(1),9000);
    x2=linspace(fitlimit(2),fitlimit(2),9000);
    y1=logspace(-5.5,4,9000);
    
    plot(x1,y1,'k'); hold on;
    plot(x2,y1,'k'); %hold on;
        

    set(gca,'XLim',[15 3e3])
    set(gca,'YLim',[7 5e3])
    set(gca,'XScale','log')
    set(gca,'YScale','log')

    xlabel('kt')
    ylabel('kt*C(kt,kt) or (kt)^3*C(kt,kt)')
    end
end
    
end