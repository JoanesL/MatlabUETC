function CohFun(pathCell,idCell,runCell,tRef,tOffSet,dx,xiscaling,tLimit)

disp(['numel ' num2str(numel(pathCell)) ])

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
if numel(tLimit)==2 && nPaths > 1
    for n = 1:nPaths
        tLim(n,:) = [tLimit(1) tLimit(2)];
    end  
else
    tLim = tLimit;
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


%pVec = 12:10:103
xlim = [0 5];
ylim = [-0.2 1.2];

for n=1:numel(pathCell)
    
    if numel(runCell{n})>1
        
        for j=1:numel(runCell{n})
            runmult=runCell{n}
        %xi Scaling turned off since ETCget does not have this option.
            [kt,r,C11ind]=UETCload(path{n},'scalar11',idCell{n},runmult(j),tRef(n),tOff{n},0);
            [kt,r,C12ind]=UETCload(path{n},'scalar12',idCell{n},runmult(j),tRef(n),tOff{n},0);
            [kt,r,C21ind]=UETCload(path{n},'scalar21',idCell{n},runmult(j),tRef(n),tOff{n},0);
            [kt,r,C22ind]=UETCload(path{n},'scalar22',idCell{n},runmult(j),tRef(n),tOff{n},0);
            [kt,r,Cvvind]=UETCload(path{n},'vector',idCell{n},runmult(j),tRef(n),tOff{n},0);
            [kt,r,Cttind]=UETCload(path{n},'tensor',idCell{n},runmult(j),tRef(n),tOff{n},0);
    
            for i=1:size(Cvvind,1)
                Cvvind(i,:) = Cvvind(i,:) .* (kt) .* (kt * r(i));
            end
        
            C11all(j,:,:)=C11ind(:,:);
            C12all(j,:,:)=C12ind(:,:);
            C21all(j,:,:)=C21ind(:,:);
            C22all(j,:,:)=C22ind(:,:);
            Cvvall(j,:,:)=Cvvind(:,:);
            Cttall(j,:,:)=Cttind(:,:);
        end
    else
        [kt,r,C11]=UETCload(path{n},'scalar11',idCell{n},runCell{n},tRef(n),tOff{n},0);
        [kt,r,C12]=UETCload(path{n},'scalar12',idCell{n},runCell{n},tRef(n),tOff{n},0);
        [kt,r,C21]=UETCload(path{n},'scalar21',idCell{n},runCell{n},tRef(n),tOff{n},0);
        [kt,r,C22]=UETCload(path{n},'scalar22',idCell{n},runCell{n},tRef(n),tOff{n},0);
        [kt,r,Cvv]=UETCload(path{n},'vector',idCell{n},runCell{n},tRef(n),tOff{n},0);
        [kt,r,Ctt]=UETCload(path{n},'tensor',idCell{n},runCell{n},tRef(n),tOff{n},0);
        for i=1:size(Cvv,1)
            Cvv(i,:) = Cvv(i,:) .* (kt) .* (kt * r(i));
        end
    
    end
        
    %R{n} = r;
    %KT{n} = kt;
%if numel(runCell{n})>1
%        UETC11{n}=C11all;
%        UETC12{n}=C12all;
%        UETC21{n}=C21all;
%        UETC22{n}=C22all;
%        UETCvv{n}=Cvvall;
%        UETCtt{n}=Cttall;
%else
%        UETC11{n}=C11;
%        UETC12{n}=C12;
%        UETC11{n}=C21;
%        UETC22{n}=C22;
%        UETCvv{n}=Cvv;
%        UETCtt{n}=Ctt;


    if numel(runCell{n})>1
        
        for j=1:numel(runCell{n})
            [k,t,E11ind]=ETCget('scalar11',idCell{n},runmult(j),tRef(n),tOff{n},0,tLim(n,:),path{n});
            [k,t,E12ind]=ETCget('scalar12',idCell{n},runmult(j),tRef(n),tOff{n},0,tLim(n,:),path{n});
            [k,t,E22ind]=ETCget('scalar22',idCell{n},runmult(j),tRef(n),tOff{n},0,tLim(n,:),path{n});
            [k,t,Evvind]=ETCget('vector',idCell{n},runmult(j),tRef(n),tOff{n},0,tLim(n,:),path{n});
            [k,t,Ettind]=ETCget('tensor',idCell{n},runmult(j),tRef(n),tOff{n},0,tLim(n,:),path{n});
            
            E11all(j,:,:)=E11ind(:,:);
            E12all(j,:,:)=E12ind(:,:);
            E22all(j,:,:)=E22ind(:,:);
            Evvall(j,:,:)=Evvind(:,:);
            Ettall(j,:,:)=Ettind(:,:);
        end
    else
            [k,t,E11]=ETCget('scalar11',idCell{n},runCell{n},tRef(n),tOff{n},0,tLim(n,:),path{n});
            [k,t,E12]=ETCget('scalar12',idCell{n},runCell{n},tRef(n),tOff{n},0,tLim(n,:),path{n});
            [k,t,E22]=ETCget('scalar22',idCell{n},runCell{n},tRef(n),tOff{n},0,tLim(n,:),path{n});
            [k,t,Evv]=ETCget('vector',idCell{n},runCell{n},tRef(n),tOff{n},0,tLim(n,:),path{n});
            [k,t,Ett]=ETCget('tensor',idCell{n},runCell{n},tRef(n),tOff{n},0,tLim(n,:),path{n});
            
    end
    
%if numel(runCell{n})>1
%        ETC11{n}=E11all;
%        ETC12{n}=E12all;
%        ETC22{n}=E22all;
%        ETCvv{n}=Evvall;
%        ETCtt{n}=Ettall;
%else
%        ETC11{n}=E11;
%        ETC12{n}=E12;
%        ETC22{n}=E22;
%        ETCvv{n}=Evv;
%        ETCtt{n}=Ett;
%end   

% * k sqrt(t t')*2 - new definition of vector


%len_kt = size(kt,2);
%len_r = length(r);


%for i=1:len_kt
%    for j=1:size(r,1)
%        R1(j,i)=r(j);
%        R2(j,i)=1/r(j);
%        Z1(j,i)=kt(i)*sqrt(r(j));
%        Z2(j,i)=kt(i)*sqrt(r(j));
%    end
%end



% for i=1:len_kt
%     for j=1:size(r,1)
%         ii = floor(i*r(j)); 
%         if ii>len_kt
%             ii = len_kt;
%         end
%         %disp(ii)
%         E11(j,i) = sqrt(C11(1,ii)*C11(1,i));
%         E12(j,i) = sqrt(C12(1,ii)*C12(1,i));
%         E22(j,i) = sqrt(C22(1,ii)*C22(1,i));
%         Ev(j,i) = sqrt(Cv(1,ii)*Cv(1,i));
%         Et(j,i) = sqrt(Ct(1,ii)*Ct(1,i));
%         
%     end
% end
% 
% E11 = ones(size(C11(:,1)))*C11(1,:);
% E12 = -ones(size(C12(:,1)))*C12(1,:);
% E22 = ones(size(C22(:,1)))*C22(1,:);
% Ev = ones(size(Cv(:,1)))*Cv(1,:);
% Et = ones(size(Ct(:,1)))*Ct(1,:);



    %DO NOT USE THIS, THE ERROR BARS IT GIVES ARE NOT THE CORRECT ONES,
    %SINCE THE DIFFERENCE BETWEEN 2 DIFFERENT RUNS OF THE SAME SIMULATION
    %IS KXI AND NOT Dii. THEREFORE THE ERROR BARS SHOULD BE CALCULATED FOR
    %KXI'S. THE EASIEST SOLUTION IS TO TAKE EVERY RUN AS AN INDIVIDUAL
    %SIMULATION.
    if numel(runCell{n})>1
        
        for j=1:numel(runCell{n})
    
            D11all(j,:,:) = C11all(j,:,:)./E11all(j,:,:);
            D12all(j,:,:) = -C12all(j,:,:)./E12all(j,:,:);
            D21all(j,:,:) = -C21all(j,:,:)./E12all(j,:,:);
            D22all(j,:,:) = C22all(j,:,:)./E22all(j,:,:);
            Dvvall(j,:,:) = Cvvall(j,:,:)./Evvall(j,:,:);
            Dttall(j,:,:) = Cttall(j,:,:)./Ettall(j,:,:);
        end
        
        D11(:,:) = mean(D11all,1);
        D12(:,:) = mean(D12all,1);
        D21(:,:) = mean(D21all,1);
        D22(:,:) = mean(D22all,1);
        Dvv(:,:) = mean(Dvvall,1);
        Dtt(:,:) = mean(Dttall,1);
        
        SDevD11=zeros(size(D11,1),size(D11,2));
        SDevD12=zeros(size(D11,1),size(D11,2));
        SDevD21=zeros(size(D11,1),size(D11,2));
        SDevD22=zeros(size(D11,1),size(D11,2));
        SDevDvv=zeros(size(D11,1),size(D11,2));
        SDevDtt=zeros(size(D11,1),size(D11,2));
        
        for j=1:numel(runCell{n})
            for k=1:size(D11,1)
                for l=1:size(D11,2)
                    SDevD11(k,l) = SDevD11(k,l) + (D11(k,l) - D11all(j,k,l)).^2;
                    SDevD12(k,l) = SDevD12(k,l) + (D12(k,l) - D12all(j,k,l)).^2;
                    SDevD21(k,l) = SDevD21(k,l) + (D21(k,l) - D21all(j,k,l)).^2;
                    SDevD22(k,l) = SDevD22(k,l) + (D22(k,l) - D22all(j,k,l)).^2;
                    SDevDvv(k,l) = SDevDvv(k,l) + (Dvv(k,l) - Dvvall(j,k,l)).^2;
                    SDevDtt(k,l) = SDevDtt(k,l) + (Dtt(k,l) - Dttall(j,k,l)).^2;
                end
            end
        end
        
        for k=1:size(D11,1)
                for l=1:size(D11,2)
                    SDevD11(k,l) = sqrt(SDevD11(k,l)/(numel(runCell{n})));
                    SDevD12(k,l) = sqrt(SDevD12(k,l)/(numel(runCell{n})));
                    SDevD21(k,l) = sqrt(SDevD21(k,l)/(numel(runCell{n}))); 
                    SDevD22(k,l) = sqrt(SDevD22(k,l)/(numel(runCell{n})));
                    SDevDvv(k,l) = sqrt(SDevDvv(k,l)/(numel(runCell{n})));
                    SDevDtt(k,l) = sqrt(SDevDtt(k,l)/(numel(runCell{n})));
                end
        end 
        
        COH11{n}=D11all;
        COH12{n}=D12all;
        COH21{n}=D21all;
        COH22{n}=D22all;
        COHvv{n}=Dvvall;
        COHtt{n}=Dttall;
        
        COH11Mean{n}=D11;
        COH12Mean{n}=D12;
        COH21Mean{n}=D21;
        COH22Mean{n}=D22;
        COHvvMean{n}=Dvv;
        COHttMean{n}=Dtt;
        
        SDevCOH11{n}=SDevD11;
        SDevCOH12{n}=SDevD12;
        SDevCOH21{n}=SDevD21;
        SDevCOH22{n}=SDevD22;
        SDevCOHvv{n}=SDevDvv;
        SDevCOHtt{n}=SDevDtt;
    else
        D11 = C11./E11;
        D12 = -C12./E12;
        D21 = -C21./E12;
        D22 = C22./E22;
        Dvv = Cvv./Evv;
        Dtt = Ctt./Ett;
    end

    k=kt/tRef(n);
    t_sim=r*tRef(n);
    
    ktmin=12.45;
    ktmax=13;
    

%HEMEN JARRI XI RUN BAKOITZARI APLIKATUZ

    if xiscaling==1              
        if numel(runCell{n})>1
            for i=1:numel(runCell{n})
                [xiLag tStat] = statGet('xiLag',idCell{n},runmult(i),path{n});
            
                which = find(t_sim<tStat(end));
                t_sim = t_sim(which);
            
                xiScale = interp1(tStat,xiLag,tRef(n));
                xiScale_sim = interp1(tStat,xiLag,t_sim);
                
                kXi(i,:) = xiScale*k;
                r_Xi(i,:) = xiScale_sim/xiScale;
            
                whichkXi{i}=find(kXi(i,:)<ktmax & kXi(i,:)>ktmin)
                kXiCell{i}=kXi(i,whichkXi{i})
                rXiCell{i}=r_Xi;
            
            %kxi(i,:)=ktind(i,:);            
            end
        %n bakoitzak ere bere kt
            KT{n}=kXiCell; 
            R{n}=rXiCell;
            disp(['//////////*********/////// ' num2str(size(KT)) ])
            disp(['//////////*********/////// ' num2str(size(KT{1})) ])

            whichKT{n}=whichkXi;
        else
            JOJOJO=size(D11)
        
            [xiLag tStat] = statGet('xiLag',idCell{n},runCell{n},path{n});
        
            which=find(t_sim<tStat(end));
            t_sim=t_sim(which);
            D11=D11(which,:);
            D12=D12(which,:);
            D21=D21(which,:);
            D22=D22(which,:);
            Dvv=Dvv(which,:);
            Dtt=Dtt(which,:);
        
            xiScale = interp1(tStat,xiLag,tRef(n));
            xiScale_sim = interp1(tStat,xiLag,t_sim);
                
            kXi=xiScale*k;
            r_Xi=xiScale_sim/xiScale;
    
            whichkXi=find(kXi<ktmax & kXi>ktmin);
            kXi=kXi(whichkXi) 
        
            KT{n}=kXi;
            R{n}=r_Xi;
            whichKT{n}=whichkXi;
            COH11{n}=D11;
            COH12{n}=D12;
            COH21{n}=D21;
            COH22{n}=D22;
            COHvv{n}=Dvv;
            COHtt{n}=Dtt;
            disp(['//////////*********/////// ' num2str(size(KT)) ])
            disp(['//////////*********/////// ' num2str(size(KT{1})) ])
        end
    else
        whichKT{n}=find( kt>ktmin & kt<=ktmax);
        kt=kt(whichkXT{n})

    end 
    %disp(['Alpha: ' alpha ])
end

%n eta HAU SOLUZIONATU BEHAR DA NOLABAIT...
%MULTIPLE RUNS BALDIN BADITTO 1.A IZAN BEHAR DA ORAIN.

%colour = {'r','m','g','c','b','b','b','b','b','b','b'};
colour = {'r','r','r','r','b','b','b','b','b','b','b'};

ax=multiPlot([2 3]);

axes(ax(1))
for n=1:numel(pathCell)
    zeinkt=whichKT{n};
    if numel(runCell{n})>1
        D11=COH11Mean{n};
        
        SDevD11=SDevCOH11{n};
    else
        D11=COH11{n};
    end
plot((R{n}-1)*KT{n},D11(:,zeinkt(1)),colour{n});
hold on;
    if numel(runCell{n})>1
    plotErrorBars((R{n}-1)*KT{n},D11(:,zeinkt(1)),SDevD11(:,zeinkt(1)),colour{n})
    hold on;
    end
end
set(gca,'XLim',[-5 5])
set(gca,'YLim',ylim)
xlabel('k(\xi(t^*)-\xi(t))')
ylabel('D_{11}') 

axes(ax(2))
for n=1:numel(pathCell)
    zeinkt=whichKT{n};
    if numel(runCell{n})>1
        D12=COH12Mean{n};
        
        SDevD12=SDevCOH12{n};
    else
        D12=COH12{n};
    end
plot((R{n}-1)*KT{n},D12(:,zeinkt(1)),colour{n});

hold on;
    if numel(runCell{n})>1
    plotErrorBars((R{n}-1)*KT{n},D12(:,zeinkt(1)),SDevD12(:,zeinkt(1)),colour{n})
    hold on;
    end
end
set(gca,'XLim',xlim)
set(gca,'YLim',ylim)
xlabel('k(\xi(t^*)-\xi(t))')
ylabel('D_{12}')

axes(ax(3))
for n=1:numel(pathCell)
    zeinkt=whichKT{n};
    
    if numel(runCell{n})>1
        D21=COH21Mean{n};
        
        SDevD21=SDevCOH21{n};
    else
        D21=COH21{n};
    end
plot((R{n}-1)*KT{n},D21(:,zeinkt(1)),colour{n});
hold on;
    if numel(runCell{n})>1
    plotErrorBars((R{n}-1)*KT{n},D21(:,zeinkt(1)),SDevD21(:,zeinkt(1)),colour{n})
    hold on;
    end
end
set(gca,'XLim',xlim)
set(gca,'YLim',ylim)
xlabel('k(\xi(t^*)-\xi(t))')
ylabel('D_{21}')
    
axes(ax(4))
for n=1:numel(pathCell)
    zeinkt=whichKT{n};
    
    if numel(runCell{n})>1
        D22=COH22Mean{n};
        
        SDevD22=SDevCOH22{n};
    else
        D22=COH22{n};
    end
plot((R{n}-1)*KT{n},D22(:,zeinkt(1)),colour{n});
hold on;
    if numel(runCell{n})>1
    plotErrorBars((R{n}-1)*KT{n},D22(:,zeinkt(1)),SDevD22(:,zeinkt(1)),colour{n})
    hold on;
    end
end
set(gca,'XLim',xlim)
set(gca,'YLim',ylim)
xlabel('k(\xi(t^*)-\xi(t))')
ylabel('D_{22}')
    
axes(ax(5))
for n=1:numel(pathCell)
    zeinkt=whichKT{n};
    
    if numel(runCell{n})>1
        Dvv=COHvvMean{n};
        
        SDevDvv=SDevCOHvv{n};
    else
        Dvv=COHvv{n};
    end
plot((R{n}-1)*KT{n},Dvv(:,zeinkt(1)),colour{n});
hold on;
    if numel(runCell{n})>1
    plotErrorBars((R{n}-1)*KT{n},Dvv(:,zeinkt(1)),SDevDvv(:,zeinkt(1)),colour{n})
    hold on;
    end 
end
set(gca,'XLim',xlim)
set(gca,'YLim',ylim)
xlabel('k(\xi(t^*)-\xi(t))')
ylabel('D_{vv}')
    
axes(ax(6))
for n=1:numel(pathCell)
    zeinkt=whichKT{n};
    
    if numel(runCell{n})>1
        Dtt=COHttMean{n};
        
        SDevDtt=SDevCOHtt{n};
    else
        Dtt=COHtt{n};
    end
plot((R{n}-1)*KT{n},Dtt(:,zeinkt(1)),colour{n});
hold on;
    if numel(runCell{n})>1
    plotErrorBars((R{n}-1)*KT{n},Dtt(:,zeinkt(1)),SDevDtt(:,zeinkt(1)),colour{n})
    hold on;
    end 
end
set(gca,'XLim',xlim)
set(gca,'YLim',ylim)
xlabel('k(\xi(t^*)-\xi(t))')
ylabel('D_{tt}')

multiPlotZoom(ax);

end

% figure
% 
% KT1 = ones(size(r)) * kt;
% KT2 = r * kt;
% 
% pcolor(KT1,KT2,log(abs(Cv)))
% shading interp
% %set(gca,'XScale','log')
% %set(gca,'YScale','log')
% set(gca,'ZScale','log')

%Delta = zeros(5,len_kt);
%dr = diff(r);
%for i=1:len_kt
%    for j=1:len_r-1
%        Delta(1,i) = Delta(1,i) + kt(i)*C11(j,i)*dr(j);
%        Delta(2,i) = Delta(2,i) + kt(i)*C12(j,i)*dr(j);
%        Delta(3,i) = Delta(3,i) + kt(i)*C22(j,i)*dr(j);
%        Delta(4,i) = Delta(4,i) + kt(i)*Cvv(j,i)*dr(j);
%        Delta(5,i) = Delta(5,i) + kt(i)*Ctt(j,i)*dr(j);
%    end
%end

%Delta = 2*Delta;

%figure(2)
% Should do label trick everywhere ..
%label = ['11' '12' '22' 'vv' 'tt'];
%for cr = 1:5
%    subplot(5,1,cr)
%    hold on
%   plot(kt,abs(Delta(cr,:)))
%    eval(['plot(kt,abs(Delta(cr,:))./abs(C' label((2*cr-1):(2*cr)) '(1,:)))'])
%    set(gca,'XScale','log')
%    set(gca,'YScale','log')
%    set(gca,'XLim',[1 1000])
%    set(gca,'YLim',[0 20])
%    hold off
%end
