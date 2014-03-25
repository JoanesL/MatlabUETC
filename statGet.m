% extract a stat from statsFile
% 2013 Version for comparing runs at different resolution
%
%Usage: stat = statGet(statName,id,run,path)
%
% path = cell array of paths to file, including final '/'
%        (if omitted gets paths from gpath global variable)

function [stat t] = statGet(statName,id,run,inPath)

if nargin==0; 
  help statGet
  return
end

global gpath

if ~exist('inPath','var'); inPath=''; end 

if numel(inPath)>0; 
  if ~ischar(inPath) 
      disp('inPath must be a string')
      return
  end
  path=inPath; 
else
  if numel(gpath)>0
  if ~ischar(gpath) 
      disp('gpath must be a string')
      return
  end
    path=gpath;
  else
    disp(['Please set gpath global variable to default path(s)'...
	  ' or specify path(s) in function call'])
    return
  end
end

%disp(['Current working path: ' path]);

%======================
%Form num2str ID string
%======================
s=strfind(id,'%');
if numel(s) > 0
    ID=[id(1:s+1) '.' id(s+1) 'i' id(s+2:end)];
else
    ID = id;
end
    
    
%==============
%Loop over runs
%==============
for r=1:numel(run)

  %=========
  %Read file
  %=========
%    if (numel(ID) > 0 && numel(run)>1) %HEMEN run=1 bada gaizki itten do... PENTSATU!!!
     file=[path 'statsFile' num2str(run(r),ID) '.dat'];
%    else
%       file=[path 'statsFile' ID '.dat'];
%    end
  disp(' ')
  %disp(['Loading file: ' file]); 
  stats=load(file);

%  disp(['File size: ' num2str(size(stats))])
  
  sizes=size(stats);
  %disp(['File size: ' num2str(sizes  )]);

  if ~exist('stat','var')
      stat = zeros(numel(run),sizes(1));
  end
  
  if sizes(2)==21
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
      eF0i_lag=stats(:,11);
      eFij=stats(:,12);
      eFij_lag=stats(:,13);
      ePi=stats(:,14);
      ePi_lag=stats(:,15);
      eDjPhi=stats(:,16);
      eDjPhi_lag=stats(:,17);
      eV=stats(:,18);
      eV_lag=stats(:,19);
      SLexp=stats(:,20);
      SLwind=stats(:,21); 
      
      
  elseif sizes(2)==16

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
  end
  
  %===============================
  %Perform Additional Calculations
  %===============================


  e=eF0i+eFij+ePi+eDjPhi+eV;    %Energy density
  L=eF0i-eFij+ePi-eDjPhi-eV;    %Lag. density
  p=eF0i-eFij/3+ePi - eDjPhi/3 - eV; % Pressure
  
  if(sizes(2) == 21)
  e_lag=eF0i_lag+eFij_lag+ePi_lag+eDjPhi_lag+eV_lag;    %Energy density
  L_lag=eF0i_lag-eFij_lag+ePi_lag-eDjPhi_lag-eV_lag;    %Lag. density
  p_lag=eF0i_lag-eFij_lag/3+ePi_lag - eDjPhi_lag/3 - eV_lag; % Pressure
  end
  
  mu = 2*pi; % mass per unit length (beta = 1 only)
  
  T0iEst_gauge = sqrt(eF0i.*eFij); % Estimate of T0i gauge
  T0iEst_scalar = sqrt(ePi.*eDjPhi); % Estimate of T0i scalar
  T0iEst = sqrt(eF0i.*eFij) + sqrt(ePi.*eDjPhi); % Estimate of T0i
  
  if(sizes(2) == 21)
  T0iEst_gauge_lag = sqrt(eF0i_lag.*eFij_lag); % Estimate of T0i gauge lag weighted
  T0iEst_scalar_lag = sqrt(ePi_lag.*eDjPhi_lag); % Estimate of T0i lag wighted
  T0iEst_lag = sqrt(eF0i_lag.*eFij_lag+ePi_lag.*eDjPhi_lag); % Estimate of T0i lag wighted
  end

switch statName
    case 'modPhi'
        stat(r,:)=stats(:,4);
    case 'modPi'
        stat(r,:)=stats(:,5);
    case 'modFij'
        stat(r,:)=stats(:,6);
    case 'modF0i'
        stat(r,:)=stats(:,7);
    case 'gaussCD'
        stat(r,:)=stats(:,8);
    case 'gaussDE'
        stat(r,:)=stats(:,9);
    case 'eF0i'
        stat(r,:)=stats(:,10);
    case 'eF0i_lag'
        stat(r,:)=stats(:,11);
    case 'eFij'
        stat(r,:)=stats(:,12);
    case 'eFij_lag'
        stat(r,:)=stats(:,13);
    case 'ePi'
        stat(r,:)=stats(:,14);
    case 'ePi_lag'
        stat(r,:)=stats(:,15);
    case 'eDjPhi'
        stat(r,:)=stats(:,16);
    case 'eDjPhi_lag'
        stat(r,:)=stats(:,17);
    case 'eV'
        stat(r,:)=stats(:,18);
    case 'eV_lag'
        stat(r,:)=stats(:,19);
    case 'SLexp'
        stat(r,:)=stats(:,20);
    case 'SLwind'
        stat(r,:)=stats(:,21); 
    case 'e'
        stat(r,:)=e; 
    case 'p'
        stat(r,:)=p; 
    case 'L'
        stat(r,:)=L; 
    case 'e_lag'
        stat(r,:)=e_lag; 
    case 'p_lag'
        stat(r,:)=p_lag; 
    case 'L_lag'
        stat(r,:)=L_lag; 
    case 'xiLag'
        stat(r,:)=sqrt(-mu./L); 
    otherwise
        error('Unknown stat')
end

end