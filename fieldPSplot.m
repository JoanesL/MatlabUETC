% field Power Spectrum plotting function for LAH data 2013
% adapted from ETCplot
%
% Usage: fieldPSplot(Fname,id,run,tOffSet,tLimit,inPath,parent)
%
%  Fname = field name, eg. E2, B2, Phi2, or dotPhi
%  id = ID string between 'fieldPS' and before '.txt'
%        eg. to load ps_B26L01 ID is '6L%2' with run=1
%  run = realizations(s) to include
%  tOffset = time when xi=0, if '*' get from statsFile (not working yet)
%  tLimit = times between which to plot
%  inPath = path to file, including final '/'
%          (if omitted or '' gets path from gpath global variable)
%  parent = destination axes for plot
%
% Version 1.1 2013.6.10 MBH

function [totalPower t]=fieldPSplot(Fname,id,run,tOffSet,tLimit,inPath,parent)

if nargin==0; 
  help fieldPSplot
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
	  ' or specify path in fucntion call'])
    return
  end
end

if ~exist('tOffSet','var'); tOffSet=0; end
if ~exist('tLimit','var'); tLimit=[0 9999999]; end

%Prepare for plot
if exist('parent','var')~=1; clf; else axes(parent); end

%Get number of runs
nRuns=size(run,2);

% %Get tOffset from statsFile if necessary
% if strcmp(tOffSet,'*')==1
%   disp(['** Getting tOffSet from statsFile Lag. fit for ' ...
% 	'tRef -> 2*tRef **'])
%   tOffSet = statsFile(-1,id,run,tRef*[1 2]);
% end
% 
% %Duplicate tOffSet if single value (eg. 0) given for many runs
% if size(tOffSet,2)==1 && nRuns>1
%   tOffSet = ones(1,nRuns)*tOffSet
% end

%Load field PS
[kk,t,Cc,sd]=fieldPSload(path,Fname,id,run,tOffSet,tLimit(1),tLimit(2));


k = kk(2:end);

nt = size(Cc,1);
nk = size(k,1);

C = zeros(nt,nk);


if ( strcmp(Fname,'B2') || strcmp(Fname,'E2') )
    for n = 1:nt
        CcRow = reshape(Cc(n,:),3,length(k));
        C(n,:) = sum(CcRow,1)';
    end
else
    C = Cc;
end
    

if exist('number','var')==1
   C=C(number(1):number(2),:);
   t=t(number(1):number(2));
end

p=1;

%Plot on top of each other (hopefully)
for i=2:size(C,1)-1 %ceil(linspace(1,size(C,1),10)) %floor(size(C,1)/2):size(C,1)
    plot(k*t(i),C(i,:)*t(i)^p,'.-g'); hold on;
end
hold on

%Highlight lines and plot uncertainties
plotLogNegative(k*t(end),C(end,:)*t(end)^p,'b',2) 
if nRuns>1 
    plotErrorBars(k*t(end),C(end,:)*t(end)^p,sd(end,:)*t(end)^p,'b')
    plotErrorBars(k*t(1),C(1,:)*t(1)^p,sd(1,:)*t(1)^p,'k')
end
plotLogNegative(k*t(1),C(1,:)*t(1)^p,'k',2);

totalPower = zeros(nt,1);
for i=1:nt
    totalPower(i) = trapz(k,k.^2.*C(i,:)'/(2*pi^2));
end


%set(gca,'XLim',[min(kt)*min(r) 40])
%set(gca,'YLim',[1e-4 1e3])
axis tight
set(gca,'XScale','log')
%set(gca,'YScale','log')
set(gca,'YScale','linear')

xlabel('kt')
ylabel('PS(kt)*t')
title(['field=' Fname])
