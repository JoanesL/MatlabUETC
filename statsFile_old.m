%StatsFile reader for LAH code 2005-2008
%
%Usage: tOffset = statsFile(pNum,id,run,tFit,dx,N,path)
%
%pNum = 0 for all or plot number for single axes set
%       eg. 2 for |phi| variation, 6 for xi
%       -1 for no plots, output Lagrangian tOffset only 
%       -2 for winding (S&V) tOffset only
%  id = ID string between 'statsFile_' and before '.dat'
%       eg. to load statsFile_6L01 ID is '6L%2' with run=1
% run = realizations(s) to include
%
%Optional parameters are:
%
% tFit = xi fit range (only needed for xi fitting)
%   dx = lattice spacing (only needed for winding xi)
%    N = lattice size (only needed for winding xi)
% path = path to file, including final '/'
%        (if ommited gets path from gpath global variable)

function tOffset = statsFile_old(pNum,id,run,tFit,dx,N,inPath)

if nargin==0; 
  help statsFile
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

disp(['Current working path: ' path])

if ~exist('dx','var'); dx=0; end
if ~exist('N','var'); N=0; end
if ~exist('tFit','var'); tFit=0; end

%==============
%Prepare figure
%==============
if pNum==0
  clf
  ax=multiPlot([3 3]);
elseif pNum>0
  clf
end

%======================
%Form num2str ID string
%======================
s=strfind(id,'%');
ID=[id(1:s+1) '.' id(s+1) 'i' id(s+2:end)];

%==============
%Loop over runs
%==============
for i=1:prod(size(run))

  %=========
  %Read file
  %=========
  file=[path 'statsFile_' num2str(run(i),ID) '.dat'];
  disp(['Loading file: ' file]); 
  stats=load(file);

  disp(['File size: ' num2str(size(stats))])

  %====================
  %Extract data columns
  %====================
  t=stats(:,1);
  a=stats(:,2);
  q=stats(:,3);
  modPhi=stats(:,4);
  modPi=stats(:,5);
  modFij=stats(:,6);
  modF0i=stats(:,7);
  gaussCD=stats(:,8);
  gaussDE=stats(:,9);
  eF0i=stats(:,10);
  eFij=stats(:,11);
  ePi=stats(:,12);
  eDjPhi=stats(:,13);
  eV=stats(:,14);
  SLexp=stats(:,15);
  SLwind=stats(:,16);

  %===============================
  %Perform Additional Calculations
  %===============================

  if prod(size(tFit))==1; tFit=[tFit max(t)]; end

  e=eF0i+eFij+ePi+eDjPhi+eV;    %Energy density
  L=eF0i-eFij+ePi-eDjPhi-eV;    %Lag. density

  mu = 2*pi;                    %Assumed string energy per unit length
  xiLag = 1./sqrt(-L/mu);

  if dx>0 & N>0
    SLlag = -(L/mu)*(N*dx)^3;
    xiWind = 1./sqrt(SLwind/(N*dx)^3);
  else
    SLlag = zeros(size(SLwind));
    xiWind = zeros(size(SLwind));
  end

  if tFit>0 & ( pNum==-1 | pNum==0 | pNum==5 | pNum==6 )
    tSub=t(t>=tFit(1) & t<=tFit(2));
    xiSub=xiLag(t>=tFit(1) & t<=tFit(2));
    pLag=leastSquares(tSub,xiSub);
    dXi_dtLag(i)=pLag(1);
    tXi0lag(i)=-pLag(2)/pLag(1);
    if prod(size(run)) == 1; 
      disp(['Lagrangian xi fitting parameters for t > ' num2str(tFit)]); 
      disp(['Gradient: ' num2str(pLag(1)) '  xi-intercept: ' num2str(pLag(2)) '  t-intercept: ' num2str(tXi0lag(i))])
    end
%    if pNum==-1
      tOffset(i) = tXi0lag(i);
%    end
  else
    pLag=[0 0];
    tXi0lag(i)=0;
  end
  tScaleLag = t - tXi0lag(i);

  if tFit>0 & dx>0 & N>0 & ( pNum==-2 | pNum==0 | pNum==5 | pNum==6 )
    pWind=leastSquares(t(t>=tFit(1) & t<=tFit(2)),xiWind(t>=tFit(1) & t<=tFit(2)));
    dXi_dtWindSV(i)=pWind(1)/sqrt((pi/6));
    tXi0wind(i)=-pWind(2)/pWind(1);
    if prod(size(run)) ==1
      disp(['Winding xi (S&V corrected) fitting parameters for t > ' num2str(tFit)])
      disp(['Gradient: ' num2str(pWind(1)/sqrt((pi/6))) '  xi-intercept: ' num2str(pWind(2)/sqrt((pi/6))) '  t-intercept: ' num2str(tXi0wind(i))])
    end
%    if pNum==-2
      tOffset(i) = tXi0wind(i);
%    end
  else
    dXi_dtWindSV(i)=0;
    tXi0wind(i)=0;
    pWind=[0 0];
  end
  tScaleWind = t - tXi0wind(i);

  %===========================
  %Output initial value of T00
  %===========================

  disp(['Initial value of T00 = ' num2str(e(1)) ' (normally 24)'])

  %============
  %PLOT RESULTS
  %============
  
  if pNum==0
    axes(ax(1))
  end
  if pNum==0 | pNum==1
    plot(t,1./(a.*q),'r'); hold on
    set(gca,'YLim',[0 max(a)])
    xlabel('time')
    ylabel('r')
    axis tight
  end

  if pNum==0
    axes(ax(2))
  end
  if pNum==0 | pNum==2
    plot(t,modPhi,'.-b',t,modPi,'.-r'); hold on
    xlabel('time')
    ylabel('|\phi| and |\pi|')
    axis tight
    hold on
  end

  if pNum==0
    axes(ax(3))
  end
  if pNum==0 | pNum==3
    plot(t,modFij./a,'c',t,modF0i./a,'m'); hold on
    xlabel('time')
    ylabel('B and E')
    axis tight
  end
    
  if pNum==0
    axes(ax(4))
  end
  if pNum==0 | pNum==4
    plot(t,gaussCD,'b',t,gaussDE,'.r',t,abs(gaussDE-gaussCD),'--k')
    set(gca,'yscale','log')
    xlabel('time')
    ylabel('Gauss Law')
    axis tight
    hold on
  end
  
  if pNum==0
    axes(ax(5))
  end
  if pNum==0 | pNum==5
     %plot(t,(xiWind-(pWind(1)*t+pWind(2)))/sqrt(pi/6),'-r'); hold on
     %plot(t,xiLag-(pLag(1)*t+pLag(2)),'-b')
     %plot([t(1) t(end)],[0 0 ],'-k')
     plot(t(2:end),diff(xiWind)./diff(t)/sqrt(pi/6),'--r'); hold on
     plot(t(2:end),diff(xiLag)./diff(t),'k'); hold on
     xlabel('\tau')
     ylabel('d\xi/d\tau')
     set(gca,'XLim',[tFit(1)*0.75 max(t)])
     plot(tFit(1)*[1 1],get(gca,'YLim'),'-k')
     plot(tFit(2)*[1 1],get(gca,'YLim'),'-k')
  end

  if pNum==0
    axes(ax(6))
  end
  if pNum==0 | pNum==6
    plot(t,xiWind,'r'); hold on;
    plot(t,xiWind/sqrt(pi/6),'--r'); %Scherrer & Vilenkin correction
    plot(t,xiLag,'b'); hold on;
    plot([0 t(end)],pWind(1)*[0 t(end)]+pWind(2),'--','Color',[0.5 0.5 0.5])
    plot([0 t(end)],pLag(1)*[0 t(end)]+pLag(2),'--','Color',[0.5 0.5 0.5])
    axis tight
    ylim=get(gca,'YLim');    
    set(gca,'YLim',[0 ylim(2)])
    xlabel('time')
    ylabel('\xi')
    plot(tFit(1)*[1 1],get(gca,'YLim'),'--','Color',[0.5 0.5 0.5])
    plot(tFit(2)*[1 1],get(gca,'YLim'),'--','Color',[0.5 0.5 0.5])
  end

  if pNum==0
    axes(ax(7))
  end
  if pNum==0 | pNum==7
    plot(t,1./sqrt(eF0i),'m',t,1./sqrt(eFij),'c',t,1./sqrt(ePi),'r',t,1./sqrt(eDjPhi),'g',t,1./sqrt(eV),'b',t,1./sqrt(e),'k')
    hold on
    xlabel('time')
    ylabel('1/sqrt(T_\mu_\nu)')
    set(gca,'YLim',[0 max(1./sqrt(0.1*e))])
  end

  if pNum==0
    axes(ax(8))
  end
  if pNum==0 | pNum==8
    plot(t,eF0i./e,'m',t,eFij./e,'c',t,ePi./e,'r',t,eDjPhi./e,'g',t,eV./e,'b')
    xlabel('time')
    ylabel('Energy fraction')
    axis tight
    hold on
  end

  if pNum==0
    axes(ax(9))
  end
  if pNum==0 | pNum==9
    plot(t,SLlag,'b',t,SLwind,'r',t,SLexp,'g'); hold on
    plot(t,SLwind*pi/6,'--r')      %Scherrer & Vilenkin corerction
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    ylabel('L / V')
    if tFit>0
      set(gca,'XLim',tFit)
    end
  end
end

if pNum==0
  multiPlotZoom(ax);
end

if prod(size(run))>1 & (pNum==-1 | pNum==0 | pNum==6)
  disp(['Mean dXi/dt = ' num2str(mean(dXi_dtLag)) ' +/- ' num2str(std(dXi_dtLag)) '  (Lag.)'])
  disp(['Mean tXi=0 = ' num2str(mean(tXi0lag)) ' +/- ' num2str(std(tXi0lag)) '  (Lag.)'])
end
if prod(size(run))>1 & (pNum==-2 | pNum==0 | pNum==6)
    disp(['Mean dXi/dt = ' num2str(mean(dXi_dtWindSV)) ' +/- ' num2str(std(dXi_dtWindSV)) '  (S&V winding.)'])
    disp(['Mean tXi=0 = ' num2str(mean(tXi0wind)) ' +/- ' num2str(std(tXi0wind)) '  (Winding.)'])
end

