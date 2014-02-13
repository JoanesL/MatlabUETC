function statsFileMulti(pNum,idCell,runCell,tFit,dx,N,pathCell,parent)
%StatsFile reader for LAH code 2005-2008
% 2013 Version for comparing runs at different resolution
%
%Usage: tOffset = statsFileMulti(pNum,idCell,runCell,tFit,dx,N,pathCell,parent)
%
%pNum     = 0 for all or plot number for single axes set
%           eg. 2 for |fields|, 6 for xi
%           -1 for no plots, output Lagrangian tOffset only 
%           -2 for winding (S&V) tOffset only
%  idCell = Cell array of ID strings between 'statsFile_' and before '.dat'
%           eg. to load statsFile_6L01 ID is '6L%2' with run=1
% runCell = Cell array of realizations(s) to include
%
%Optional parameters are:
%
% tFit     = xi fit range (only needed for xi fitting)
%   dx     = lattice spacing (only needed for winding xi)
%    N     = lattice size (only needed for winding xi)
% pathCell = cell array of paths to file, including final '/'
%        (if omitted gets paths from gpath global variable)
%
% Version 1.1 2013.7.25 MBH

if nargin==0; 
  help statsFileMulti
  return
end

global gpath

if ~exist('pathCell','var'); pathCell={}; end 

if numel(pathCell)>0; 
  if ~iscell(pathCell) 
      disp('pathCell must be a cell array')
      return
  end
  path=pathCell; 
else
  if numel(gpath)>0
  if ~iscell(gpath) 
      disp('gpath must be a cell array')
      return
  end
    path=gpath;
  else
    disp(['Please set gpath global variable to default path(s)'...
	  ' or specify path(s) in function call'])
    return
  end
end

if ~exist('dx','var'); dx=0; end
if ~exist('N','var'); N=0; end
if ~exist('tFit','var'); tFit=0; end

%Duplicate parameters if single value given for many paths
nPaths = numel(path);
if numel(idCell)==1 && nPaths > 1
    for n = 1:nPaths
        idCell(n) = idCell(1);
    end  
end
if numel(runCell)==1 && nPaths > 1
    for n = 1:nPaths
        runCell(n) = runCell(1);
    end  
end
if numel(dx)==1 && nPaths > 1
    for n = 1:nPaths
        dx(n) = dx(1);
    end  
end
if numel(N)==1 && nPaths > 1
    for n = 1:nPaths
        N(n) = N(1);
    end  
end
if numel(tFit)==1 && nPaths > 1
    for n = 1:nPaths
        tFit(n) = tFit(1);
    end  
end

%==============
%Prepare figure
%==============
if pNum==0
    if ~exist('parent','var');
      clf
      ax=multiPlot([3 3]);
    else
      ax=parent;
    end
elseif pNum>0
    if ~exist('parent','var'); 
        clf; 
        ax=gca
    else
        ax=axes(parent); 
    end
end

%==============
%Plot figure(s)
%==============

for n=1:numel(path)

    id = idCell{n};
    run = runCell{n};
    inPath = path{n};

    statsFile(pNum,id,run,tFit(n),dx(n),N(n),inPath,ax);

end

if pNum==0
    multiPlotZoom(ax);
end

