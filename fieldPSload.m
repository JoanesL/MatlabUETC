%function [k,t,Fav,sd]=fieldPSload(filenamePre,Fname,id,run,tRef,tOffSet,tStart,tEnd)
function [k,t,Cav,sd]=fieldPSload(filenamePre,Fname,id,run,tOffSet,tStart,tEnd)

%Get number of runs to load
nRuns=size(run,2);

%Make id into num2str format
s=strfind(id,'%');
if numel(s) > 0
    ID=[id(1:s+1) '.' id(s+1) 'i' id(s+2:end)];
    idString = num2str(run(1),ID);
else
    ID = id;
    idString = ID;
end

%Load info file for first number
[k,tSim]=fieldPSinfoRead(filenamePre,idString);

% %Store simulation times
% tSim=tRef*r;
% tSimRef=tRef;
% k=kt/tSimRef;

tSimRef = tSim(1);

%Load correlators
for i=1:nRuns
      if numel(s) > 0
      	idString = num2str(run(i),ID);
      else
        idString = ID;
      end
    fileName=[filenamePre 'ps_' Fname idString '.txt'];
    disp(['Loading field PS: ' fileName]);
    Ctemp=load(fileName);
    if size(tSim,1)<size(Ctemp,1); 
        disp('Warning: correlator has more rows than expected. Truncating...')
        Ctemp=Ctemp(1:size(tSim,1),:); 
    end
    if size(tSim,1)>size(Ctemp,1); 
        disp('Warning: correlator has less rows than expected...')
        disp('Removing extra times in tSim vector')
        disp(['Initial size of t = ' num2str(size(r,1))])
        tSim=tSim(1:size(Ctemp,1));
        C=C(1:size(Ctemp,1),:,:);
        disp(['New size of tSim = ' num2str(size(tSim,1))])
    end
    if i==1
        C=Ctemp;
    else	
        C(:,:,i)=Ctemp;
    end

    %Apply time correction
    if i==1
        T=tSim-tOffSet(i);
    else
        T(:,i)=tSim-tOffSet(i);
    end
    for j=1:size(T,1);
    	C(j,:,i)=C(j,:,i)*tSim(j)/T(j,i);
    end
end 

%Interpolate onto grid
if std(tOffSet)>0
   disp('Interpolating runs onto grid')
   t=linspace(max(T(1,:)),min(T(end,:)),size(r,1));
   
   [Ki,Ti]=meshgrid(k,t);
   
   for i=1:nRuns
       Ci(:,:,i)=interp2(k,T(:,i),C(:,:,i),Ki,Ti,'cubic');
   end
else
   t=T(:,1);
   Ci=C;
end   

%Average power spectra
for i=1:nRuns
   if i==1; 
        Cav=Ci(:,:,1);
        Csqu=Ci(:,:,1).^2;
    else
        Cav=Cav+Ci(:,:,i);
        Csqu=Csqu+Ci(:,:,i).^2;
    end
end

Cav=Cav/nRuns;
if nargout==4
    if nRuns>1
        sd=sqrt(Csqu/nRuns-Cav.^2)*sqrt(nRuns/(nRuns-1));   %SD from sample
        sd=sd/sqrt(nRuns);              %SD on mean
    else
       sd=zeros(size(Cav));
   end
end

%Remove early and late times
which=t+mean(tOffSet)>=tStart & t+mean(tOffSet)<=tEnd;
t=t(which);
Cav=Cav(which,:);
if nargout==4
  sd=sd(which,:);
end
disp(['Min/max simulation times: ' num2str(min(t+mean(tOffSet))) '  ' num2str(max(t+mean(tOffSet)))])

