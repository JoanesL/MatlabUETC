function [k,t,C,sd]=ETCplot(Cname,id,run,tRef,tOffSet,tLimit,inPath,parent)
%ETC plotting function for UETC.hpp data 2006-2008
% revisions 2013 MBH 
% - plotting without scaling (useful for standing wave sims)
% - plotting with xi scaling
%
%Usage: ETCplot(Cname,id,run,tRef,tOffset,tLimit,path,parent)
%
%Cname = ETC name, eg. scalar11 or vector
%   id = ID string between 'ETCscalar11_' and before '.dat'
%        eg. to load statsFile_6L01 ID is '6L%2' with run=1
%  run = realizations(s) to include
% tRef = UETC reference time
% tOffset = time when xi=0, 
% - if '*' get from statsFile (Lag & tRef<t<2*tRef)
% - if 'noscaling' do not scale k with time (and undo correlator time scaling) 
% - if 'xiscaling' get xi from statsFile using statGet
%
%Optional parameters are:
%
%  tLimit = vector with times between which to plot
%   path = path to file, including final '/'
%          (if omited or '' gets path from gpath global variable)
% parent = destination axes for plot

if nargin==0; 
  help ETCplot
  return
end

global gpath

if ~exist('inPath','var'); inPath=''; end 

if numel(inPath)>0; 
  path=inPath; 
else
  if numel(gpath)>0
    path=gpath;
  else
    disp(['Please set gpath global variable to default path'...
	  ' or specify path in function call'])
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
  tOffSet = statsFile(-1,id,run,tRef*[1 2],0.5,1024,path); % dx, N kludge
end

%Check to see if we want to scale
noscaling = 0;
if strcmp(tOffSet,'noscaling')==1
  disp('** Plotting against k and removing t factor from correlator')
  tOffSet = 0;
  noscaling = 1;
end

%Get xiLag from statsFile 
xiscaling = 0;
if strcmp(tOffSet,'xiscaling')==1
  disp(['** Scaling with xiLag'])
  [xiLag tStat] = statGet('xiLag',id,run,path);
  if nRuns > 1
      xiLagAv = mean(xiLag,1);
  else
      xiLagAv = xiLag;
  end
  tOffSet = 0;
  xiscaling = 1;
end


%Duplicate tOffSet if single value (eg. 0) given for many runs
if size(tOffSet,2)==1 && nRuns>1
  tOffSet = ones(1,nRuns)*tOffSet;
end


%Load ETC

[k,t,C,sd]=ETCload(path,Cname,id,run,tRef,tOffSet,tLimit(1),tLimit(2));
C=abs(C);

mint = min(t);
maxt = max(t);

if exist('number','var')==1
   C=C(number(1):number(2),:);
   t=t(number(1):number(2));
end

if (noscaling == 1)
    for i=1:size(C,1)
        C(i,:) = C(i,:)/t(i);
        t(i) = 1;
    end
end

if (xiscaling == 1)
    %This must be done since usually tStat(end)<t(end), therefore it is
    %imposible to perform the interpolation to the lastest times.
    which=find(t<tStat(end));
    t=t(which);
    C=C(which,:);
    for i=1:size(C,1)
        xiScale = interp1(tStat,xiLagAv,t(i));
        if strcmp(Cname,'vector')~=1
            C(i,:) = xiScale*C(i,:)/t(i);
        else
            C(i,:) = C(i,:)*(t(i)/xiScale);
        end
        t(i) = xiScale;
    end
end


%Plot on top of each other (hopefully)
for i=2:size(C,1)-1 %ceil(linspace(1,size(C,1),10)) %floor(size(C,1)/2):size(C,1)
    plot((k*t(i)),C(i,:),'.-g'); hold on;
end
hold on

t(end-2)

%Highlight lines and plot uncertainties
plotLogNegative((k*t(end)),C(end,:),'b',2) 
if nRuns>1 
    plotErrorBars((k*t(end)),C(end,:),sd(end,:),'b')
    plotErrorBars((k*t(1)),C(1,:),sd(1,:),'k')
end
plotLogNegative((k*t(1)),C(1,:),'k',2);

%set(gca,'XLim',[min(kt)*min(r) 40])
%set(gca,'YLim',[1e-4 1e3])
axis tight
set(gca,'XScale','log')
%set(gca,'YScale','log')
set(gca,'YScale','linear')

if (noscaling==1)
    xlabel('k')
    ylabel('C(kt,kt)/t')
elseif (xiscaling==1)
    xlabel('k\xi')
    ylabel('C(k\xi,k\xi)/\xi')
else    
    xlabel('kt')
    ylabel('C(kt,kt)')
end
title(['C=' Cname ', min time' num2str(mint) ', max time' num2str(maxt)])
