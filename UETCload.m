function [kt,r,C]=UETCload(filenamePre,Cname,id,run,tRef,tOffSet,xiscaling)
%UETC loading function for UETC.hpp data 2006-2013
%
% Usage: [kt,r,C]=UETCload(filenamePre,Cname,id,run,tRef,tOffset)
%
% filenamePre = path to data
%  Cname = UETC name, eg. scalar11 or vector
%     id = ID string between 'UETCscalar11_' and before '.dat'
%          eg. to load statsFile_6L01 ID is '6L%2' with run=1
%    run = realizations(s) to include
%   tRef = UETC reference time
%tOffset = time when xi=0, if '*' get from statsFile (Lag & tRef<t<2*tRef)
%xiscaling = 1-> xi scaling, toffset=0
%            0-> toffset rescaling
%
% 13.10.04 MBH tOffset now performed for each run rather than on average
%              Can now handle data files without runID
%
% TO DO: allow xi scaling (as with ETCplot/load) in UETCtimeOffset


%Get number of runs to load
nRuns=size(run,2);

if xiscaling==1
    tOffSet=0;
end

%Make id into num2str format
s=strfind(id,'%');
for nr=1:nRuns
    if s > 0
        ID=[id(1:s+1) '.' id(s+1) 'i' id(s+2:end)];
        idString{nr} = num2str(run(nr),ID);
    else
        ID = id;
        idString{nr} = ID;
    end
end


if exist('tOffSet','var')
    %Duplicate tOffSet if single value (eg. 0) given for many runs
    if(xiscaling==0)
        if strcmp(tOffSet,'*')==1
            tOffSet=tOffSet;
        elseif size(tOffSet,2)==1 && nRuns>1
            tOffSet = ones(1,nRuns)*tOffSet;
        end        
        
    end
end

%Load info file for first number 
[ktRaw,rRaw]=UETCinfoRead(filenamePre,idString{1});

%Load and sum correlators
for i=1:nRuns
    fileName=[filenamePre 'UETC' Cname idString{i} '.dat'];
    disp(['Loading UETC: ' fileName]);
    Ci=load(fileName);
    
    if size(Ci,1)>size(rRaw,1); 
        disp('Warning correlator has more rows than expected. Truncating...')
        Ci=Ci(1:size(rRaw,1),:);
    end
    
    if numel(tOffSet)>1 && nRuns>1
        timeOffSet = tOffSet(i);
    else
        timeOffSet = tOffSet;
    end

%    if exist('tOffSet','var')==1
        disp(['i ' num2str(i) ' toffset ' num2str(timeOffSet) ])
        disp(['Initial max time ratio: ' num2str(max(rRaw))])
        [kt,r,Ci]=UETCtimeOffSet(Cname,id,run,ktRaw,rRaw,Ci,tRef,timeOffSet,xiscaling,filenamePre); %Apply time offset correction
                
        disp(['Post-offset max time ratio: ' num2str(max(r))])
    
    if i==1; 
        C=Ci;
    else
        C=C+Ci;
    end
end

%Average
if nRuns>1
    disp('Averaging correlators');
    C=C/nRuns;
end

