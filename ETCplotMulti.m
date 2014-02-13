function ETCplotMulti(correl,idCell,runCell,tRef,tOffset,tLimit,pathCell,scale)
%
% ETCplotMulti(correl,id,run,tRef,tOffset,tLimit,pathCell,scale)
%
% Plots ETCs from several directories, on the same graph.
%
%Cname = ETC name, [scalar11|scalar12|scalar22|vector|tensor]
%  idCell = Cell array of ID strings between 'statsFile_' and before '.dat'
%           eg. to load statsFile_6L01 ID is '6L%2' with run=1
% runCell = Cell array of realizations(s) to include
%
%    tRef = UETC reference time
% tOffset = time when xi=0, 
% - if '*' get from statsFile (Lag & tRef<t<2*tRef)
% - if 'noscaling' do not scale k with time (and undo correlator time scaling) 
% - if 'xiscaling' get xi from statsFile using statGet
%
%  tLimit = vector with times between which to plot
%  pathCell = cell array of paths to files, including final '/'
% e.g. pathCell = {'1024/uetc/','2048/uetc/','4096/uetc/'};
%  scale = 'log' (default) or 'linear'
%
% example:
% 
% p = {'prod_LAH_stat/Mat_all/uetc/','s_tref_test/mat_s10/uetc/'};
% ETCplotMulti('scalar11',{'_mat_%2',''},{5:7,1},300,'xiscaling',[300 330],p)
%
% Version 1.1 2013.6.11 MBH
% Version 1.2 2013.7.30 MBH

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
if numel(tRef)==1 && nPaths > 1
    for n = 1:nPaths
        tRef(n) = tRef(1);
    end  
end
if numel(tLimit)==2 && nPaths > 1
    for n = 1:nPaths
        tLim(n,:) = [tLimit(1) tLimit(2)];
    end  
else
    tLim = tLimit;
end
% Convert to cell array even for one element
for n = 1:nPaths
    tOff{n} = tOffset;
end  

if ~exist('scale','var')
    scale = 'log';
end

clf

for n=1:numel(pathCell)
    disp(tOff(n))
    h = subplot(1,1,1);
    ETCplot(correl,idCell{n},runCell{n},tRef(n),tOff{n},tLim(n,:),path{n},h);
    set(h,'YScale',scale)

end
