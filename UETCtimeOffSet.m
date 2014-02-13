function [kt,r,C]=UETCtimeOffSet(Cname,id,run,kt,r,C,tRef,tOffSet,xiscaling,path)
%UETC scaling function for UETC.hpp data 2006-2013
%
% Usage: [kt,r,C]=UETCtimeOffset(Cname,kt,r,C,tRef,tOffset)
%
%  Cname = UETC name, eg. scalar11 or vector
%     kt = wavenumber*time
%      r = time'/time
%      C = Correlator array
%   tRef = UETC reference time
%tOffset = time when xi=0, if '*' get from statsFile (Lag & tRef<t<2*tRef)
%xiscaling = 1-> xi scaling, toffset=0
%            0-> toffset rescaling
%
% 13.10.04 MBH tOffset now performed for each run rather than on average
%              Can now handle data files without runID
%
% TO DO: allow xi scaling (as with ETCplot/load) in UETCtimeOffset
%        to be implemented when strcomp(tOffset,'xiscaling') is true

nRuns=size(run,2);

%Get tOffset from statsFile if necessary

if strcmp(tOffSet,'*')==1
  disp(['** Getting tOffSet from statsFile Lag. fit for ' ...
    'tRef -> 2*tRef **'])
  tOffSet = statsFile(-1,id,run,[tRef (tRef*(4/3))],0.5,4096,path); % Kludge to get path in
  tOffSet=mean(tOffSet);
end

%Get xiLag from statsFile 
if xiscaling==1
  disp(['** Scaling with xiLag'])

  [xiLag tStat] = statGet('xiLag',id,run,path);
    
  if nRuns > 1
      xiLagAv = mean(xiLag,1);
  else
      xiLagAv = xiLag;
  end
  tOffSet = 0;
end

%Adjust times
tSim=tRef*r;
tSimRef=tRef;
t=tSim-tOffSet;
tRef=tRef-tOffSet;

%Form outputs
r=t./tRef;
kt=kt.*tRef/tSimRef;

if (xiscaling == 1)
    %for i=1:size(C,1)
        xiScale = interp1(tStat,xiLagAv,t);
        %t(i) = xiScale;
    %end

end

if (xiscaling == 1)
    %Scale UETC with Xi(tRef) and Xi(t')
    if strcmp(Cname,'vector')~=1
        C=C.*repmat(sqrt((xiScale./t))*sqrt((xiScale(1)/tRef)),size(kt));
    else
        C=C.*repmat(sqrt((t./xiScale))*sqrt((tRef/xiScale(1))),size(kt));
    end
else
    if nargout>2
        if strcmp(Cname,'vector')~=1
        C=C.*repmat(sqrt(t./tSim)*sqrt(tRef/tSimRef),size(kt));
        else
        C=C.*repmat(sqrt(tSim./t)*sqrt(tSimRef/tRef),size(kt));
        disp('** This version has correct vector tOffset corretion **')
        end
    end
end



