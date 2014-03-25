function UETCmerge(pathCell,outputpath,idCell,runCell,era,tRef,tOffSet,dx,xiscaling,highkXiCorrection,NewNormalization,WidthNormalization,LowkXiMode)

%function UETCmerge(Cname,pathCell,idCell,runCell,tRef,tOffSet,dx,xiscaling,normalization,prog,NewNormalization,WidthNormalization)

FirstUETC=1;
NumUETC=1;

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
if numel(dx)==1 && nPaths > 1
    for n = 1:nPaths
        dx(n) = dx(1);
    end  
end
% Convert to cell array even for one element
for n = 1:nPaths
    tOff{n} = tOffSet;
end  

if ~exist('scale','var')
    scale = 'log';
end

%Get number of runs
%nRuns=size(run,2);

%Get tOffset from statsFile if necessary
if strcmp(tOffSet,'*')==1
    disp(['** Getting tOffSet from statsFile Lag. fit for ' ...
	'tRef -> 2*tRef **'])
    for n=1:numel(pathCell)
        tOffSet(n) = statsFile(-1,idCell{n},runCell{n},tRef(n),dx(n),4096,path{n});
    end
end


%Duplicate tOffSet if single value (eg. 0) given for many runs
%if size(tOffSet,2)==1 && nRuns>1
%  tOffSet = ones(1,nRuns)*tOffSet
%end
for n=1:numel(pathCell)
    if n==1
        disp('#####################')
        disp('S=0')
        disp('#####################')
    elseif n==2
        disp('#####################')
        disp('S=1')
        disp('#####################')
    end
    for i=FirstUETC:NumUETC
        if i==1
            Cname='scalar11';
        elseif i==2
            Cname='scalar22';
        elseif i==3
            Cname='vector';
        elseif i==4
            Cname='tensor';
        elseif i==5
            Cname='scalar12';
        end
        disp('@@@@@@@@@@@@@')
        disp('Loading UETCs..')
        disp('@@@@@@@@@@@@@')
        [kt,r,C1]=UETCload(path{n},Cname,idCell{n},runCell{n},tRef(n),tOff{n},xiscaling);
      
        %figure()
        %plot(kt,abs(C1(1,:)),'b')
        %set(gca,'XScale','log')
        %ylabel('JODER')
        %xlabel('k(\xi \xi`)^{1/2}')
   
        if highkXiCorrection==1
            disp('@@@@@@@@@@@@@')
            disp('Loading ETCs..')
            disp('@@@@@@@@@@@@@')
            [k,t,E11]=ETCget(Cname,idCell{n},runCell{n},tRef(n),tOff{n},1,[0 9000],path{n});
        SIZEE11=size(E11)
        end
          
        k=kt/tRef(n);
        t_sim=r*tRef(n);
   
        if strcmp(Cname,'scalar12')==1
            [kt,r,C12]=UETCload(path{n},'scalar21',idCell{n},runCell{n},tRef(n),tOff{n},xiscaling);
        end
   %if strcmp(Cname,'vector')==1
   %     for i=1:size(C1,1)
   %         C1(i,:) = C1(i,:) .* (kt) .* (kt * r(i)); 
   %     end
   %end
      
        disp('@@@@@@@@@@@@@')
        disp('Getting kXi and r_xi...')
        disp('@@@@@@@@@@@@@')
        
        if(xiscaling==1) %Behin eginda nahikoa
            if numel(runCell{n})>1
                runmult=runCell{n};
                for l=1:numel(runCell{n})
                    [xiLag(l,:) tStat] = statGet('xiLag',idCell{n},runmult(l),path{n});
                    xiLagS=xiLag(l,:);
                end
        
                %Remove late times due to velocity decay
                if n==2
                    if strcmp(era,'mat')==1
                        which=find(t_sim<800);
                    elseif strcmp(era,'rad')==1
                        which=find(t_sim<600);
                    end
                else
                    which=find(t_sim<tStat(end));
                end
                t_sim=t_sim(which);
        
                xiLagAv=mean(xiLagS,1);        
        
                xiScale = interp1(tStat,xiLagAv,tRef(n));
                xiScale_sim = interp1(tStat,xiLagAv,t_sim);
                
                kXi=xiScale*k;
                r_Xi=xiScale_sim/xiScale;
            else
                [xiLagAv tStat] = statGet('xiLag',idCell{n},runCell{n},path{n});
        
                if n==2
                    which=find(t_sim<800);
                else
                    which=find(t_sim<tStat(end));
                end
        
                t_sim=t_sim(which);
                
                C1=C1(which,:);
                if strcmp(Cname,'scalar12')==1
                    C1_12=C1_12(which,:);
                end
        
                xiScale = interp1(tStat,xiLagAv,tRef(n));
                xiScale_sim = interp1(tStat,xiLagAv,t_sim);
                
                kXi=xiScale*k;
                r_Xi=xiScale_sim/xiScale;
            end
            
            %Behin jakindakoan, korrekzioa aplikatu
            C1=C1(which,:);
            if strcmp(Cname,'scalar12')==1
                C12=C12(which,:);
            end
            if highkXiCorrection==1
                E11=E11(which,:);
            end
            
        else
            r_Xi = t_sim/tRef(n);
            kXi = kt;
        end
        
        if n==1
            UETCS0{i}=C1;
            if strcmp(Cname,'scalar12')==1
                UETCS0{i+1}=C12;
            end
            if highkXiCorrection==1
                ETCS0{i}=E11;
            end
        elseif n==2
            UETCS1{i}=C1;
            if strcmp(Cname,'scalar12')==1
                UETCS1{i+1}=C12;
            end
            if highkXiCorrection==1
                ETCS1{i}=E11;
            end
        end
    end
    %s=0 eta s=1-ei dagokiona
    R{n} = r_Xi;
    KXi{n} = kXi;
end

%Equal for all correlators
    r1=R{1};
    r2=R{2};
    
    kt1=KXi{1};
    kt2=KXi{2};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DEFINITION OF VARIABLES AND PARAMETERS
    
    %for interpolation
    rlimit=r2(end)
    
    rmerges1_Size=ceil(150*((rlimit-1)/(r1(end)-1)))
    rmerges0_Size=150-rmerges1_Size;
    
    rmerges1=linspace(1,rlimit,rmerges1_Size);
    
    resMerge=(r1(end)-1)/150;
    rmerges0=linspace(rlimit+resMerge,r1(end),rmerges0_Size);
    
    rmerge=[rmerges1 rmerges0];
    
    %MAIN LOOP
for i=FirstUETC:NumUETC
    if i==1
        Cname='scalar11';
    elseif i==2
        Cname='scalar22';
    elseif i==3
        Cname='vector';
    elseif i==4
        Cname='tensor';
    elseif i==5
        Cname='scalar12';
    end
        
    %LOAD CORRELATORS
    C1=UETCS0{i};
    C2=UETCS1{i}; 
    
    UETCS0{i}=[];
    
    if strcmp(Cname,'scalar12')==1
        C1_12=UETCS0{i+1};
        C2_12=UETCS1{i+1};
        UETCS0{i+1}=[];
    end
       
    %Localize UETC peak in k\xi's
          
    if strcmp(era,'rad')==1
        if strcmp(Cname,'scalar11')==1
            PeakkXi=find(kt1<6 & kt1>1.75);
        elseif strcmp(Cname,'scalar12')==1
            PeakkXi=find(kt1<8 & kt1>3);  
        elseif strcmp(Cname,'scalar22')==1
            PeakkXi=find(kt1<7 & kt1>1);
        elseif strcmp(Cname,'vector')==1
            PeakkXi=find(kt1<3 & kt1>1);
        elseif strcmp(Cname,'tensor')==1
            PeakkXi=find(kt1<10 & kt1>1);
        end
    elseif strcmp(era,'mat')==1
        if strcmp(Cname,'scalar11')==1
            PeakkXi=find(kt1<6 & kt1>1.5);
        elseif strcmp(Cname,'scalar12')==1
            PeakkXi=find(kt1<7 & kt1>3);  
        elseif strcmp(Cname,'scalar22')==1
            PeakkXi=find(kt1<8 & kt1>1);
        elseif strcmp(Cname,'vector')==1
            PeakkXi=find(kt1<3 & kt1>1);
        elseif strcmp(Cname,'tensor')==1
            PeakkXi=find(kt1<10 & kt1>1);
        end
    end
    
    Peak_kXi{i}=PeakkXi;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %///////////////////////////////////
    %HIGH kXi Correction. 
    %///////////////////////////////////
    
    %It should be done first since uses ETCs at original r values.
    
    if highkXiCorrection==1
        ETC11_1=ETCS0{i};
        if strcmp(Cname,'scalar12')==1
            ETC1=ETCS0{1};
            ETC2=ETCS0{2};
        C1=UETCextrap(C1,kt1,r1,ETC1,ETC2);
        C1_12=UETCextrap(C1_12,kt1,r1,ETC1,ETC2);
        else
        C1=UETCextrap(C1,kt1,r1,ETC11_1);
        end
              
        ETC11_2=ETCS1{i}; 
        if strcmp(Cname,'scalar12')==1
            ETC1=ETCS1{1};
            ETC2=ETCS1{2};
        C2=UETCextrap(C2,kt1,r1,ETC1,ETC2);
        C2_12=UETCextrap(C2_12,kt1,r1,ETC1,ETC2);
        else
        C2=UETCextrap(C2,kt2,r2,ETC11_2);
        end
    end
        
    
    %ADD zeros AT LOW kt
    %LATER EXTRAPOLATE DATA
    if min(kt1)<min(kt2)
        kt2extra=min(kt1):0.5*min(kt1):kt2(1);
        kt2_extra=[kt2extra kt2(1:end)];
        
        diffkt2=size(kt2extra);
        
        Zero=zeros(size(C2,1));
        
        C2=[repmat(Zero(:,1),size(kt2extra)) C2(:,1:end)];
                
        if strcmp(Cname,'scalar12')==1
            C2_12=[repmat(C2_12(:,1),size(kt2extra)) C2_12(:,1:end)];
        end
                  
    end
    
    %INTERPOLATE UETC TO SELECTED KT's    
    for l=1:size(C2,1)
        C2i(l,:)=interp1(kt2_extra,C2(l,:),kt1);
    end
    if strcmp(Cname,'scalar12')==1
        for l=1:size(C2_12,1)
            C2i_12(l,:)=interp1(kt2_extra,C2_12(l,:),kt1); 
        end
    end
    
    UETCS1_i{i}=C2i;
    if strcmp(Cname,'scalar12')==1
        UETCS1_i{i+1}=C2i_12;
    end
        
    UETCS0_Orig{i}=C1;
    if strcmp(Cname,'scalar12')==1
        UETCS0_Orig{i+1}=C1_12;
    end
    
    %///////////////////////////////////
    %ETC AND WIDTH NORMALIZATION
    %///////////////////////////////////
    
    if NewNormalization==1
        %Perform ETC normalization, for values higher than extra kt
        %normETC=C2i(1,size(which_extrakXi,2)+1:end)./C1(1,size(which_extrakXi,2)+1:end);

        %Calculate factor only for s=0 kt values
        normETC=C2i(1,3+1:end)./C1(1,3+1:end);
        
        FactorETC{i}=normETC;
        
        figure()
        ax1=multiPlot([1 2]);
        axes(ax1(1))
        
        plot(kt1,abs(C1(1,:)),'b',kt1,abs(C2i(1,:)),'r')
        set(gca,'XScale','log')
        ylabel('(\xi/t)|C^s_{11}(kt,1)|')
        xlabel('k(\xi \xi`)^{1/2}')      
        
        %Normalize all s=0 UETC using this correction
        for l=1:size(C1,1)
            %C1(i,size(which_extrakXi,2)+1:end)=normETC.*C1(i,size(which_extrakXi,2)+1:end);
            C1(l,6+1:end)=normETC((3+1):end).*C1(l,6+1:end);
        end
            if strcmp(Cname,'scalar12')==1
                for l=1:size(C1_12,1)
                    C1_12(l,6+1:end)=normETC((3+1):end).*C1_12(l,6+1:end);
                end
            end
            
        UETCS0{i}=C1;
        if strcmp(Cname,'scalar12')==1
            UETCS0{i+1}=C1_12;
        end    
        axes(ax1(2))
        
        plot(kt1(3+1:end),normETC,'b')
        set(gca,'XScale','log')
        xlabel('k(\xi \xi`)^{1/2}')
        ylabel('\gamma (k\xi)')
        
        %Bi hok kendu!!!
    end
end

%         for l=1:size(C2i,2)
%             C2Norm(:,l)=interp1(r2,C2i(:,l),rmerges1);
%         end
%         UETCS1_Norm{i}=C2Norm;
%             if strcmp(Cname,'scalar12')==1
%                 for l=1:size(C2i_12,2)
%                     C2Norm_12(:,l)=interp1(r2,C2i_12(:,l),rmerges1);
%                 end
%                 UETCS1_Norm{i+1}=C2Norm_12;
%             end
%         
%         for l=1:size(C1,2)
%             C1NormTemp(:,l)=interp1(r1,C1(:,l),rmerge);
%         end    
%             if strcmp(Cname,'scalar12')==1
%                 for l=1:size(C1_12,2)
%                     C1NormTemp_12(:,l)=interp1(r1,C1_12(:,l),rmerge);
%                 end
%             end
%         %%%
%     
%         %Something we'll need later to add low kXi data
%         %first 6 points normalized with the mean value of ETCnorm
%         C1_LowkXi=C1NormTemp(:,1:(6));
%         UETC_LowkXi{i}=C1_LowkXi;
%         if strcmp(Cname,'scalar12')==1
%             C1_LowkXi_12=C1NormTemp_12(:,1:(6));
%             UETC_LowkXi{i+1}=C1_LowkXi_12;
%         end
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%
%             
%         %///////////////////////////////////
%         %WIDTH NORMALIZATION
%         %///////////////////////////////////
%         
%           if WidthNormalization==1
%           
%           C1=UETCS0{i};
%           if strcmp(Cname,'scalar12')==1
%               C1_12=UETCS0{i+1};
%           end
%           
%           C2Norm=UETCS1_Norm{i};
%           if strcmp(Cname,'scalar12')==1
%               C2Norm_12=UETCS1_Norm{i+1};
%           end
%           
%           PeakkXi=Peak_kXi{i};
%               
% %           if i==1
% %               figure()
% %               ax2=multiPlot([2 3]);
% %           end
%           
%           axes(ax2(i))
%           %for j=1:size(PeakkXi,2)
%           for j=10:10
%             plot(r1,abs(C1(:,PeakkXi(j))),'b',rmerges1,abs(C2Norm(:,PeakkXi(j))),'r')
%             hold on;
%           end
%             set(gca,'XScale','log')
%             xlabel('r_{\xi}')
%             ylabel('((\xi/t)(\xi`/t`))^{1/2} |C^s_{11}(ktMax)|')
%             set(gca,'XLim',[1 1.5],'XTick',[0.125 0.25 0.5 1 2 4 8])
%             %set(gca,'YLim',[0 50])
%           
%             %Find to which r corresponds the last C2 value
%           
%             rnorm=linspace(rlimit-0.1,rlimit+0.1,200);
%             
%             for l=1:size(C1,2)
%                  C1rNorm(:,l)=interp1(r1,C1(:,l),rnorm);
%             end
%                 if strcmp(Cname,'scalar12')==1
%                     for l=1:size(C1_12,2)
%                         C1rNorm_12(:,l)=interp1(r1,C1_12(:,l),rnorm);
%                     end
%                 end
%                 
%             for j=1:size(PeakkXi,2)
%                 zein=find(C1rNorm(:,PeakkXi(j))<=(C2Norm(end,PeakkXi(j))+0.2) & C1rNorm(:,PeakkXi(j))>=(C2Norm(end,PeakkXi(j))-0.2));
%                 
%                 if size(zein,1)>1
%                     rfin(j)=mean(rnorm(zein));
%                 else 
%                     rfin(j)=rnorm(zein);
%                 end
%                        
%                 %interpolate functions taking into account rfin
%                 
%                 res=(r1(end)-1)/150; %space between different r values.
%                 
%                 rs1(j,:)=linspace(1,rfin(j)-res,rmerges1_Size);
%                 rs0(j,:)=linspace(rfin(j),r1(end),rmerges0_Size);
%                                                 
%                 rWidth(j,:)=[rs1(j,:) rs0(j,:)];
%                 
%                 C1WidthNorm(:,j)=interp1(r1,C1(:,PeakkXi(j)),rWidth(j,:));
% 
%                 C1WidthPrueba(:,j)=interp1(rWidth(j,:),C1WidthNorm(:,j),rmerge);
%                 
%                 %Correct r-values greater than rfin using rFactor
%                 rFactor(j)=rfin(j)-rlimit;
%                 rs0Corrected(j,:)=rs0(j,:)-rFactor(j);
%                 
%                 %rFactor(j)=rlimit/rfin(j);
%                 %rs0Corrected(j,:)=rs0(j,:)*rFactor(j);
%                                 
%                 rWidthCorrected(j,:)=[rs1(j,:) rs0Corrected(j,:)];
%                 
%                 %Once we have the corrected r's, use those to interpolate
%                 %to a single set of r's.
% 
%                 C1WidthCorrected(:,j)=interp1(rWidthCorrected(j,:),C1WidthNorm(:,j),rmerge);
%                 
%                 %Insert r-corrected values in the original UETC matrix
%                 C1NormTemp(:,PeakkXi(j)) = C1WidthCorrected(:,j);
%                 
%                 if j==size(PeakkXi,2)
%                     rMax(i)=min(rs0Corrected(:,end))
%                     
%                     %Erase parameters
%                     rfin=[];
%                     rs1(j,:)=[];
%                     rs0(j,:)=[];
%                     rWidth=[];;
%                     C1WidthNorm=[];
%                     C1WidthPrueba=[];
%                     rFactor=[];
%                     rs0Corrected=[];
%                     rWidthCorrected=[];
%                     C1WidthCorrected=[];
%                 end
%                 
%                     if strcmp(Cname,'scalar12')==1
%                         zein_12=find(C1rNorm_12(:,PeakkXi(j))<=(C2Norm_12(end,PeakkXi(j))+0.1) & C1rNorm_12(:,PeakkXi(j))>=(C2Norm_12(end,PeakkXi(j))-0.1));
%                 
%                         if size(zein_12,1)>1
%                             rfin_12(j)=mean(rnorm(zein_12));
%                         else 
%                             rfin_12(j)=rnorm(zein_12);
%                         end
%                        
%                         %interpolate functions taking into account rfin
%                         rs1_12(j,:)=linspace(1,rfin_12(j)-0.01,size(rmerges1,2));
%                         rs0_12(j,:)=linspace(rfin_12(j),r1(end),size(rmerges0,2));
%                                                 
%                         rWidth_12(j,:)=[rs1_12(j,:) rs0_12(j,:)];
%                 
%                         C1WidthNorm_12(:,j)=interp1(r1,C1_12(:,PeakkXi(j)),rWidth_12(j,:));
% 
%                         C1WidthPrueba_12(:,j)=interp1(rWidth_12(j,:),C1WidthNorm_12(:,j),rmerge);
%                 
%                         %Correct r-values greater than rfin using rFactor
%                         rFactor_12(j)=rfin_12(j)-rlimit;
%                         rs0Corrected_12(j,:)=rs0_12(j,:)-rFactor_12(j);
%                 
%                         %rFactor_12(j)=rlimit/rfin_12(j);
%                         %rs0Corrected_12(j,:)=rs0_12(j,:)*rFactor_12(j);
%                                 
%                         rWidthCorrected_12(j,:)=[rs1_12(j,:) rs0Corrected_12(j,:)];
%                 
%                         %Once we have the corrected r's, use those to interpolate
%                         %to a single set of r's.
% 
%                         C1WidthCorrected_12(:,j)=interp1(rWidthCorrected_12(j,:),C1WidthNorm_12(:,j),rmerge);
%                 
%                         %Insert r-corrected values in the original UETC matrix
%                         C1NormTemp_12(:,PeakkXi(j)) = C1WidthCorrected_12(:,j);
%                         if j==size(PeakkXi,2)
%                             rMax_12=min(rs0Corrected_12(:,end));
%                             rMax(i+1)=rMax_12;
%                             
%                             %Erase parameters
%                             rfin_12=[];
%                             rs1_12(j,:)=[];
%                             rs0_12(j,:)=[];
%                             rWidth_12=[];;
%                             C1WidthNorm_12=[];
%                             C1WidthPrueba_12=[];
%                             rFactor_12=[];
%                             rs0Corrected_12=[];
%                             rWidthCorrected_12=[];
%                             C1WidthCorrected_12=[];
%                         end
%                     end
%             end 
%             UETCS0_Norm{i}=C1NormTemp;
%             if strcmp(Cname,'scalar12')==1
%                 UETCS0_Norm{i+1}=C1NormTemp_12;
%                 C1NormTemp_12=[];
%             end 
%             C1NormTemp=[];
%           end
%     end
%     
% end
% 
% if WidthNormalization==1
%     %Determine rMax
%     rMaxFinal=min(rMax)
%     
%     rmergeCorrected=linspace(1,rMaxFinal,150);
%     
%     for i=FirstUETC:NumUETC
%         if i==1
%             Cname='scalar11';
%         elseif i==2
%             Cname='scalar22';
%         elseif i==3
%             Cname='vector';
%         elseif i==4
%             Cname='tensor';
%         elseif i==5
%             Cname='scalar12';
%         end
%         
%         C1NormTemp=UETCS0_Norm{i};
%         UETCS0_Norm{i}=[];
%         if strcmp(Cname,'scalar12')==1
%             C1NormTemp_12=UETCS0_Norm{i+1};
%             UETCS0_Norm{i+1}=[];
%         end
%         
%         C2Norm=UETCS1_Norm{i};
%         if strcmp(Cname,'scalar12')==1
%             C2Norm_12=UETCS1_Norm{i+1};
%         end
%         
%         for j=1:size(C1NormTemp,2) 
%             C1NormTot(:,j)=interp1(rmerge,C1NormTemp(:,j),rmergeCorrected);
%         end   
%         UETCS0_Norm{i}=C1NormTot;
%         if strcmp(Cname,'scalar12')==1
%             for j=1:size(C1NormTemp_12,2) 
%                 C1NormTot_12(:,j)=interp1(rmerge,C1NormTemp_12(:,j),rmergeCorrected);
%             end
%             UETCS0_Norm{i+1}=C1NormTot_12;
%         end
%         
%         if i==1
%             figure()
%             ax3=multiPlot([2 3]);
%         end
%         
%         C1=UETCS0{i};
%         
%         PeakkXi=Peak_kXi{i};
%         axes(ax3(i))
%         for j=10:10
%             plot(rmergeCorrected,abs(C1NormTot(:,PeakkXi(j))),'b',rmerges1,abs(C2Norm(:,PeakkXi(j))),'r')
%             hold on;
%         end
%         set(gca,'XScale','log')
%         xlabel('r_{\xi}')
%         ylabel('((\xi/t)(\xi`/t`))^{1/2} |C^s_{11}(ktMax)|')
%         set(gca,'XLim',[1 1.5],'XTick',[0.125 0.25 0.5 1 2 4 8])
%     end
%     
% else
%     rmergeCorrected=rmerge;
% end
% 
%     
%     
% %///////////////////////////////////
% %MERGING
% %///////////////////////////////////
% 
%     
%     whichrOrig_s0=find(r1>=rlimit);
%     whichrCorrected_s0=find(rmergeCorrected>=rlimit);
%     whichrCorrected_s1=find(rmergeCorrected<=rlimit);
%     whichrOrig_s1=find(r2<=rlimit);
%     
%     r1s=r1(whichrOrig_s0);
%     r1sCorr=rmergeCorrected(whichrCorrected_s0);
%     r2s=r2(whichrOrig_s1);
% 
%     if WidthNormalization==1
%         rs=[r2s;r1sCorr'];
%     else
%         rs=[r2s;r1s];
%     end
%     
%     %%PRUEBA
%     for i=FirstUETC:NumUETC
%         if i==1
%             Cname='scalar11';
%         elseif i==2
%             Cname='scalar22';
%         elseif i==3
%             Cname='vector';
%         elseif i==4
%             Cname='tensor';
%         elseif i==5
%             Cname='scalar12';
%         end
%         
%         if WidthNormalization==1
%             C1=UETCS0_Norm{i};
%         else
%             C1=UETCS0{i};
%         end
%         C1Orig=UETCS0_Orig{i};
%         C2i=UETCS1_i{i};
%         C2Orig=UETCS1{i};
%             if strcmp(Cname,'scalar12')==1
%                 if WidthNormalization==1
%                     C1_12=UETCS0_Norm{i+1};
%                 else
%                     C1_12=UETCS0{i+1};
%                 end
%                 C1Orig_12=UETCS0_Orig{i+1};
%                 C2_12=UETCS1_i{i+1};
%                 C2Orig_12=UETCS1{i+1};
%             end
%             
%             
%             
%         C2prob=C2Orig(whichrOrig_s1,:);
%         C1prob=C1Orig(whichrOrig_s0,:);
%         if strcmp(Cname,'scalar12')==1
%             C2_12prob=C2Orig_12(whichrOrig_s1,:);
%             C1_12prob=C1Orig_12(whichrOrig_s0,:);
%         end
%     
%         if WidthNormalization==1
%             C1=C1(whichrCorrected_s0,:);
%         else
%             C1=C1(whichrOrig_s0,:);
%         end
%         C2=C2i(whichrOrig_s1,:);
%         if strcmp(Cname,'scalar12')==1
%             if WidthNormalization==1
%                 C1_12=C1_12(whichrCorrected_s0,:);
%             else
%                 C1_12=C1_12(whichrOrig_s0,:);
%             end
%             C2_12=C2i_12(whichrOrig_s1,:);
%         end
%     
%     %MERGE DATA
%         C2rsize=size(C2,1);
%         C1rsize=size(C1,1);
%         for l=1:C1rsize
%             C2(C2rsize+l,:)=C1(l,:);
%         end
%     
%         if strcmp(Cname,'scalar12')==1
%             C2_12rsize=size(C2_12,1);
%             C1_12rsize=size(C1_12,1);
%             for l=1:C1_12rsize;
%                 C2_12(C2_12rsize+l,:)=C1_12(l,:);
%             end
%         end
%     
%         %FINALLY INTERPOLATE TO SELECTED RATIOS
%         for l=1:size(C2,2)
%             CTOT(:,l)=interp1(rs,C2(:,l),rmergeCorrected);
%         end
%         if strcmp(Cname,'scalar12')==1
%             for l=1:size(C2_12,2)
%                 CTOT_12(:,l)=interp1(rs,C2_12(:,l),rmergeCorrected);
%             end
%         end
% 
%         rMax=rmergeCorrected(end)
%         
%     %ADD LOW KXI DATA, s=0 data 
%     
%     %Get mean value of the normalization using first 5 points
%     
%         C1_LowkXi=UETC_LowkXi{i};
%         for l=1:size(C1_LowkXi,2)
%             C1_AddLowkXi(:,l)=interp1(rmerge,C1_LowkXi(:,l),rmergeCorrected);
%             if strcmp(Cname,'scalar12')==1
%                 C1_LowkXi_12=UETC_LowkXi{i+1};
%                 C1_AddLowkXi_12(:,l)=interp1(rmerge,C1_LowkXi_12(:,l),rmergeCorrected);
%             end
%         end
%     
%         if LowkXiMode==1
%             normETC=FactorETC{i};
%             normETCmean=0;
%             l2=0;
%             for l=1:5
%                 normETCmean = normETCmean + (l*l)*normETC(l);
%                 l2 = l2 + (l*l);
%             end
%     
%             normETCmean = normETCmean/l2
%     
%             CTOT(:,1:(6))=normETCmean*C1_AddLowkXi;
%             if strcmp(Cname,'scalar12')==1
%                 CTOT_12(:,1:(6))=normETCmean*C1_AddLowkXi_12;
%             end
%         elseif LowkXiMode==2
%             CTOT(:,1:(6))=C1_AddLowkXi;
%             if strcmp(Cname,'scalar12')==1
%                 CTOT_12(:,1:(6))=C1_AddLowkXi_12;
%             end
%         end
%     %///////////////////////////////////
%     %BACK TO kt and non-Xiscaling UETC. Then OUTPUT
%     %///////////////////////////////////
%     %Use statGet and linear fit to get the mean slope of xi
%     %Always s=1 -> n=2
%         
%         alpha = statsFile(20,idCell{2},runCell{2},[tRef(2) (tRef(2)*(4/3))],0.5,4096,pathCell{2})
%         
%         %Corrections
%         ktFinal=kt1/alpha;
%         
%         if strcmp(Cname,'vector')~=1
%             CFinal=(1/alpha)*CTOT;
%         else
%             CFinal=alpha*CTOT; %Correct vector normalization
%         end
%         
%         if strcmp(Cname,'scalar12')==1
%             CFinal_12=(1/alpha)*CTOT_12;
%         end
% %             CTOT=[];
% %             if strcmp(Cname,'scalar12')==1
% %                 CTOT_12=[];
% %             end
%         
%     %OUTPUT
%     
% %     %Info files,
% %     sizekt=size(CFinal,2);
% %     sizer=size(CFinal,1);
% %     if strcmp(Cname,'scalar11')==1 %We need just one info file
% %         if strcmp(era,'mat')==1
% %             infoFilename = [outputpath  'info_mat_01.dat'];
% %         elseif strcmp(era,'rad')==1
% %             infoFilename = [outputpath  'info_rad_01.dat'];
% %         end
% %         
% %         infoCell{1}=sizekt;
% %         infoCell{2}=ktFinal;
% %         infoCell{3}=sizer;
% %         infoCell{4}=rmergeCorrected;
% %                 
% %         dlmwrite(infoFilename,'Number','delimiter','')
% %         dlmwrite(infoFilename,sizekt,'-append','precision','%i')
% %         dlmwrite(infoFilename,'ktFinal','-append','delimiter','')
% %         dlmwrite(infoFilename,ktFinal,'-append','delimiter',' ','precision','%f')
% %         dlmwrite(infoFilename,'Number','-append','delimiter','')
% %         dlmwrite(infoFilename,sizer,'-append','precision','%i')
% %         dlmwrite(infoFilename,'r','-append')
% %         dlmwrite(infoFilename,rmergeCorrected,'-append','delimiter',' ','precision','%f')
% %         
% %         %fid = fopen(infoFilename, 'w');
% % 
% %         %for row=1:4
% %         %    if row==1
% %         %        fprintf(fid, '%d\n', mycell{row,:});
% %         %    elseif row==2
% %         %        fprintf(fid, '%d\n', mycell{row,:});
% %         %end
% % 
% %         %fclose(fid);
% %         
% %         %save(infoFilename,'sizekt','ktFinal','sizer','rmergeCorrected','-ascii','-tabs');
% %     end
% %     %UETCs
% %     if strcmp(era,'mat')==1
% %         uetcFilename = [outputpath 'UETC' Cname '_mat_01.dat'];
% %         if strcmp(Cname,'scalar12')==1
% %             uetcFilename_12 = [outputpath 'UETCscalar21_mat_01.dat'];
% %             %save(uetcFilename_12,'CFinal_12','-ascii');
% %             dlmwrite(uetcFilename_12,CFinal_12,' ');
% %         end
% %     elseif strcmp(era,'rad')==1
% %         uetcFilename = [outputpath 'UETC' Cname '_rad_01.dat'];
% %         if strcmp(Cname,'scalar12')==1
% %             uetcFilename_12 = [outputpath 'UETCscalar21_rad_01.dat'];
% %             %save(uetcFilename_12,'CTOT_12','-ascii');
% %             dlmwrite(uetcFilename_12,CFinal_12,' ');
% %         end
% %     end
% %     %save(uetcFilename,'CTOT','-ascii');
% %     dlmwrite(uetcFilename,CFinal,' ');
%     
%     %///////////////////////////////////
%     %PLOT
%     %///////////////////////////////////
%     
% for i=1:size(kt1,2)
%     for j=1:size(rmergeCorrected,2)
%         R1(j,i)=rmergeCorrected(j);
%         R2(j,i)=1/rmergeCorrected(j);
%         Z(j,i)=kt1(i)*sqrt(rmergeCorrected(j));
%         Z2(j,i)=kt1(i);
%     end
% end 
% 
% %%PROBATAKO
% 
% colourprob0=colourMap(C1prob);
% colourprob1=colourMap(C2prob);
%     if strcmp(Cname,'scalar12')==1
%         colourprob0_12=colourMap(C1_12prob);
%         colourprob1_12=colourMap(C2_12prob);
%     end
% for i=1:size(kt1,2)
%     for j=1:size(r1s,1)
%         R1s0(j,i)=r1s(j);
%         R2s0(j,i)=1/r1s(j);
%         Zs0(j,i)=kt1(i)*sqrt(r1s(j));
%         Z2s0(j,i)=kt1(i);
%     end
% end 
% for i=1:size(kt2,2)
%     for j=1:size(r2s,1)
%         R1s1(j,i)=r2s(j);
%         R2s1(j,i)=1/r2s(j);
%         Zs1(j,i)=kt2(i)*sqrt(r2s(j));
%         Z2s1(j,i)=kt2(i);
%     end
% end 
% 
% figure()
% ax=multiPlot([2 1]);
% 
% axes(ax(1))
% surf(R1s0,Z2s0,abs(C1prob),colourprob0,'FaceLighting','phong','FaceColor','interp','EdgeColor','flat')
%    hold on;
% if strcmp(Cname,'scalar12')==1
% surf(R2s0,Z2s0,abs(C1_12prob),colourprob0_12,'FaceLighting','phong','FaceColor','interp','EdgeColor','flat')
% else
% surf(R2s0,Z2s0,abs(C1prob),colourprob0,'FaceLighting','phong','FaceColor','interp','EdgeColor','flat')
% end
% 
% surf(R1s1,Z2s1,abs(C2prob),colourprob1,'FaceLighting','phong','FaceColor','interp','EdgeColor','flat')
%    hold on;
% if strcmp(Cname,'scalar12')==1
% surf(R2s1,Z2s1,abs(C2_12prob),colourprob1_12,'FaceLighting','phong','FaceColor','interp','EdgeColor','flat')
% else
% surf(R2s1,Z2s1,abs(C2prob),colourprob1,'FaceLighting','phong','FaceColor','interp','EdgeColor','flat')
% end
% 
% %projection
% %view(0,0)
% view([-141 36])
% camlight right
% camlight left
% camlight headlight
% set(gca,'YScale','log','XScale','log')
% set(gca,'ZScale','linear')
%     
% axis tight
% set(gca,'YLim',[1 200])
% %set(gca,'ZLim',[10^(-3) 1])
% set(gca,'box','on')
% set(gca,'FontSize',14)
% 
% set(gcf,'Color',[1 1 1])
% set(gca,'Color',[0.8 0.8 0.8]);
% 
% set(gca,'LineWidth',2);
% set(gca,'FontSize',14)
% 
% if(xiscaling == 1)
%    xlabel('r_{\xi}')
%    ylabel('k (\xi \xi`)^{1/2}')
%    if strcmp(Cname,'scalar11')==1
%       zlabel('((\xi/t)(\xi`/t`))^{1/2} |C^s_{11}|')
%   elseif strcmp(Cname,'scalar12')==1
%        zlabel('((\xi/t)(\xi`/t`))^{1/2} |C^s_{12}|')    
%    elseif strcmp(Cname,'scalar21')==1
%        zlabel('((\xi/t)(\xi`/t`))^{1/2} |C^s_{21}|')
%    elseif strcmp(Cname,'scalar22')==1
%        zlabel('((\xi/t)(\xi`/t`))^{1/2} |C^s_{22}|')
%    elseif strcmp(Cname,'vector')==1
%        zlabel('((\xi/t)(\xi`/t`))^{1/2} |C^s_{vv}|')
%    elseif strcmp(Cname,'tensor')==1
%        zlabel('((\xi/t)(\xi`/t`))^{1/2} |C^s_{tt}|')
%    end
% else
%    xlabel('r')
%    ylabel('K (t t`)^{1/2}')
%    if strcmp(Cname,'scalar11')==1
%        zlabel('|C^s_{11}|')
%    elseif strcmp(Cname,'scalar12')==1
%        zlabel('|C^s_{12}|')   
%    elseif strcmp(Cname,'scalar21')==1
%        zlabel('|C^s_{21}|')
%    elseif strcmp(Cname,'scalar22')==1
%        zlabel('|C^s_{22}|')
%    elseif strcmp(Cname,'vector')==1
%        zlabel('|C^s_{vv}|')
%    elseif strcmp(Cname,'tensor')==1
%       zlabel('|C^s_{tt}|')
%    end
% end
% 
% rMax=8;
% rMax=rmerge(end);
% ktMin=min(Z);
% ktMin=min(ktMin);
% ktMax=max(Z,[],1);
% ktMax=max(ktMax);    
% CMax=max(CTOT);
% CMax=max(CMax);
% set(gca,'XLim',[1/rMax rMax],'XTick',[0.125 0.25 0.5 1 2 4 8])
% set(gca,'YLim',[ktMin ktMax])
% %set(gca,'ZLim',[0 CMax])
% %%%%%%%%%%%%%%%%%%%
%   
% axes(ax(2))
% %colour=colourMap(CFinal);
% colour=colourMap(CTOT);
% 
% % surf(R1,Z,abs(CFinal),colour,'FaceLighting','phong','FaceColor','interp','EdgeColor','flat')
% %     hold on;
% % if strcmp(Cname,'scalar12')==1
% %     colour2=colourMap(CFinal_12);
% %     surf(R2,Z,abs(CFinal_12),colour2,'FaceLighting','phong','FaceColor','interp','EdgeColor','flat')
% % else
% %     surf(R2,Z,abs(CFinal),colour,'FaceLighting','phong','FaceColor','interp','EdgeColor','flat')
% % end
% 
% surf(R1,Z,abs(CTOT),colour,'FaceLighting','phong','FaceColor','interp','EdgeColor','flat')
%     hold on;
% if strcmp(Cname,'scalar12')==1
%     colour2=colourMap(CTOT_12);
%     surf(R2,Z,abs(CTOT_12),colour2,'FaceLighting','phong','FaceColor','interp','EdgeColor','flat')
% else
%     surf(R2,Z,abs(CTOT),colour,'FaceLighting','phong','FaceColor','interp','EdgeColor','flat')
% end
% %projection
% %view(0,0)
% view([-141 36])
% %camlight right
% %camlight left
% camlight headlight
% set(gca,'YScale','log','XScale','log')
% set(gca,'ZScale','linear')
%     
% axis tight
% %set(gca,'YLim',[1 200])
% %set(gca,'ZLim',[10^(-3) 1])
% set(gca,'box','on')
% set(gca,'FontSize',14)
% 
% 
% if(xiscaling == 1)
%     xlabel('r_{\xi}')
%     ylabel('K (\xi \xi`)^{1/2}')
%     if strcmp(Cname,'scalar11')==1
%         zlabel('((\xi/t)(\xi`/t`))^{1/2} |C^s_{11}|')
%     elseif strcmp(Cname,'scalar12')==1
%         zlabel('((\xi/t)(\xi`/t`))^{1/2} |C^s_{12}|')    
%     elseif strcmp(Cname,'scalar21')==1
%         zlabel('((\xi/t)(\xi`/t`))^{1/2} |C^s_{21}|')
%     elseif strcmp(Cname,'scalar22')==1
%         zlabel('((\xi/t)(\xi`/t`))^{1/2} |C^s_{22}|')
%     elseif strcmp(Cname,'vector')==1
%         zlabel('((\xi/t)(\xi`/t`))^{1/2} |C^s_{vv}|')
%     elseif strcmp(Cname,'tensor')==1
%         zlabel('((\xi/t)(\xi`/t`))^{1/2} |C^s_{tt}|')
%     end
% else
%     xlabel('r')
%     ylabel('K (t t`)^{1/2}')
%     if strcmp(Cname,'scalar11')==1
%         zlabel('|C^s_{11}|')
%     elseif strcmp(Cname,'scalar12')==1
%         zlabel('|C^s_{12}|')   
%     elseif strcmp(Cname,'scalar21')==1
%         zlabel('|C^s_{21}|')
%     elseif strcmp(Cname,'scalar22')==1
%         zlabel('|C^s_{22}|')
%     elseif strcmp(Cname,'vector')==1
%         zlabel('|C^s_{vv}|')
%     elseif strcmp(Cname,'tensor')==1
%         zlabel('|C^s_{tt}|')
%     end
% end
% 
% set(gcf,'Color',[1 1 1])
% set(gca,'Color',[0.8 0.8 0.8]);
% 
% set(gca,'LineWidth',2);
% set(gca,'FontSize',14)
% 
% %rMax=8;
% rMax=rmerge(end);
% ktMin=min(Z);
% ktMin=min(ktMin);
% ktMax=max(Z,[],1);
% ktMax=max(ktMax);    
% CMax=max(CTOT);
% CMax=max(CMax);
% set(gca,'XLim',[1/rMax rMax],'XTick',[0.125 0.25 0.5 1 2 4 8])
% set(gca,'YLim',[ktMin ktMax])
% 
% %set(gca,'ZLim',[0 CMax])
%     CFinal=[];
%     if strcmp(Cname,'scalar12')==1
%         CFinal_12=[];
%     end
%     end
%     multiPlotZoom(ax)
 end

function C=colourMap(C)

amp=max(max(abs(C)));
R=zeros(size(C)); G=R; B=R;

%Create smooth red to yellow transition for positive values
bool=C>0;
RtoG=(C/amp)*(pi/4)+pi/4;
R(bool)=sin(RtoG(bool));
G(bool)=cos(RtoG(bool));

%Create green blue transition for range 1-amp down to 0
bool=(bool==0);
GtoB=-(C/amp)*(pi/2);
G(bool)=cos(GtoB(bool));
B(bool)=sin(GtoB(bool));

%Create brightest possible map
C(:,:,1)=R; C(:,:,2)=G; C(:,:,3)=B;
maxRGB=max(C,[],3);
C(:,:,1)=C(:,:,1)./maxRGB;
C(:,:,2)=C(:,:,2)./maxRGB;

end

function D=UETCextrap(C,kt,r,ETC,ETC2)

%Hard-wired setting
ktExtrap = 30; %Measured at reference time, kXi in our case.

disp(['Performing power-law correction from k\xi = ' num2str(ktExtrap) ' (hardcoded)'])
disp('Requires ETCload to use ksqrt(tt'')C rather than C')

%Find kt first index (at reference time) to correct and adjust reference kt to match
[temp,iExtrap] = min(abs(kt - ktExtrap));
ktExtrap = kt(iExtrap);

%Determine correction/Extrapolation matrix using ETC data for symmetric cases
%(assumes 1/kt dependence)
if nargin==4
  R=ones(size(C));
  R(1,iExtrap:end) = ( ETC(1,iExtrap) * ktExtrap ./ kt(iExtrap:end) ) ./ ETC(1,iExtrap:end);
  for j=2:max(size(ETC,1))
    R(j,iExtrap:end) = sqrt( R(1,iExtrap:end) .* ( ETC(j,iExtrap:end) * ktExtrap ./ kt(iExtrap:end) ) ./ ETC(j,iExtrap:end) );
  end
end

%Determine correction/extrapolation matrix using ETC and ETC2 data from unsymmetric cases
%(convention for ETC(j) and ETC2(1) matches that uses by UETC.hpp such that C12 involves)
%(phi (ie. 1) from a general time and psi (ie. 2) from tRef)
if nargin==5
  R=ones(size(C));
  for j=1:max(size(ETC,1))
    R(j,iExtrap:end) = sqrt( ( ETC(j,iExtrap) * ktExtrap ./ kt(iExtrap:end) ) ./ ETC(j,iExtrap:end) ) ...
                        .* sqrt( ( ETC2(1,iExtrap) * ktExtrap ./ kt(iExtrap:end) ) ./ ETC2(1,iExtrap:end) );
  end
end

%Apply correction/extrapolation
D=C.*R;
end


    %PERFORM NORMALIZATION OF THE S=0 SIDE AT THE BOUNDARY
    %if normalization==1
    %    s0limit=abs(C1(1,:));
    %    s1limit=abs(C2(end,:));
    %    
    %    normfactor=s1limit./s0limit;
    %    
    %    normnon=find(kt1<8);
    %    
    %    figure()
    %    ax1=multiPlot([1 2]);
    %
    %    axes(ax1(1))
    %    plot(kt1(normnon)*sqrt(r1s(1)),abs(C1(1,normnon)),'b')
    %    set(gca,'XScale','log')
    %    xlabel('k(\xi \xi`)^{1/2}')
    %    ylabel('((\xi/t)(\xi`/t`))^{1/2} |C^s_{11}(r_{\xi}^{limit})|')
    %    axes(ax1(2))
    %    plot(kt1(normnon)*sqrt(r1s(1)),normfactor(normnon),'b')
    %    set(gca,'XScale','log')
    %    xlabel('k(\xi \xi`)^{1/2}')
    %    ylabel('\gamma (k\xi)')
    %    
    %    num=40;
    %    
    %    if strcmp(Cname,'scalar12')==1
    %        s0limit2=abs(C1_12(1,:));
    %        s1limit2=abs(C2_12(end,:));
    %        normfactor2=s1limit2./s0limit2;
    %        
    %        if prog==1
    %            for i=1:num
    %                C1_12(i,normnon)=(((normfactor2(normnon)-1)/num)*(num-i+1)+1).*C1_12(i,normnon);
    %            end
    %        else
    %            for i=1:size(C1_12,1)
    %                C1_12(i,normnon)=normfactor2(normnon).*C1_12(i,normnon);
    %            end
    %        end
    %        
    %    end
    %    
    %    if prog==1
    %        for i=1:num
    %            C1(i,normnon)=(((normfactor(normnon)-1)/num)*(num-i+1)+1).*C1(i,normnon);
    %        end
    %    else
    %        for i=1:size(C1,1)
    %            C1(i,normnon)=normfactor(normnon).*C1(i,normnon);
    %        end
    %    end
    %   
    %end