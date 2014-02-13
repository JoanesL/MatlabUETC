%ETC get function for UETC.hpp data 2006-2008
%
%Usage: ETCget(Cname,id,run,tRef,tOffset,tLimit,path,parent)
%
%Cname = ETC name, eg. scalar11 or vector
%   id = ID string between 'ETCscalar11_' and before '.dat'
%        eg. to load statsFile_6L01 ID is '6L%2' with run=1
%  run = realizations(s) to include
% tRef = UETC reference time
%tOffset = time when xi=0, if '*' get from statsFile (Lag & tRef<t<2*tRef)
%
%Optional parameters are:
%
%  after = time after which to plot
%   path = path to file, including final '/'
%          (if ommited or '' gets path from gpath global variable)

function [k,t,C,sd]=ETCget(Cname,id,run,tRef,tOffSet,xiscaling,tLimit,inPath)

disp(['path ETCget-en ' inPath ])

if nargin==0; 
  help ETCget
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
if exist('parent','var')~=1; clf; else axes(parent); end

%Get number of runs
nRuns=size(run,2);

%Get tOffset from statsFile if necessary
if strcmp(tOffSet,'*')==1
  disp(['** Getting tOffSet from statsFile Lag. fit for ' ...
	'tRef -> 2*tRef **'])
  tOffSet = statsFile(-1,id,run,tRef*[1 2]);
end

%Duplicate tOffSet if single value (eg. 0) given for many runs
if size(tOffSet,2)==1 && nRuns>1
  tOffSet = ones(1,nRuns)*tOffSet
end

%Load ETC
[k,t,C,sd]=ETCload(path,Cname,id,run,tRef,tOffSet,tLimit(1),tLimit(2));
C=abs(C);
if (xiscaling == 1)
    [xiLag tStat] = statGet('xiLag',id,run,path);
    if nRuns > 1
      xiLagAv = mean(xiLag,1);
    else
      xiLagAv = xiLag;
    end
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
    end
end

