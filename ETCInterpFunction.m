function ETCInterpFunction(Cname,s,tEq,idRM,runRM,pathCell,idCell,runs,tref,tOffset,tLimit_rm,tLimit,Fit)
 
eraB=[1.09 1.17 1.44 1.76 1.91];

    for i=1:size(tEq,2)
        
        nRunsRM=size(runRM,2);
        nRuns=size(runs,2);
        
        id=['_tEq' num2str(tEq(i)) idRM];
        idIns=['_eraB' num2str(eraB(i)) idRM];
        
            %Get tOffset from statsFile if necessary
        if strcmp(tOffset,'*')==1
            disp(['** Getting tOffSet from statsFile Lag. fit for ' ...
            'tRef -> 2*tRef **'])
            tOffSet = statsFile(-1,id,runRM,150*[1 2],0.5,1024,pathCell{1}); % dx, N kludge
            tEq(i)=tEq(i)-tOffSet;
        end

        %Check to see if we want to rescale
        noscaling = 0;
        if strcmp(tOffset,'noscaling')==1
            disp('** Plotting against k and removing t factor from correlator')
            tOffSet_temp=0;
            tOffSet = ones(1,nRunsRM)*tOffSet_temp;
            noscaling = 1;
        end
        
        %Get xiLag from statsFile 
        xiscaling = 0;
        if strcmp(tOffset,'xiscaling')==1
            disp(['** Scaling with xiLag'])
            [xiLag tStat] = statGet('xiLag',id,runRM,pathCell{1});
            A = statsFile(20,id,runRM,150*[1 2],0.5,1024,pathCell{1}); % dx, N kludge
            B = statsFile(-1,id,runRM,150*[1 2],0.5,1024,pathCell{1});
            xiEq(i)=A*(tEq(i)-B);

                xiLagAv = xiLag;
                
                tOffSet = 0;
                xiscaling = 1;
                tEq(i)=xiEq(i);
        end
        
        %Load Correlators
        [k_rm,t,C_rm,sd_rm]=ETCload(pathCell{1},Cname,id,runRM,150,tOffSet,tLimit_rm(1),tLimit_rm(2));
        %[kInst,tInst,C_rmInst,sd_rm]=ETCload(pathCell{1},Cname,idIns,runRM,150,tOffSet,150,150.5);
         
        %t_teqInst(i)=tInst/tEq(i);
        
        if (xiscaling == 1)
            %This must be done since usually tStat(end)<t(end), therefore it is
            %imposible to perform the interpolation to the lastest times.
            which=find(t<tStat(end));
            t=t(which);
            C_rm=C_rm(which,:);
            for l=1:size(C_rm,1)
                xiScale = interp1(tStat,xiLagAv,t(l));
                if strcmp(Cname,'vector')~=1
                    C_rm(l,:) = xiScale*C_rm(l,:)/t(l);
                else
                    C_rm(l,:) = C_rm(l,:)*(t(l)/xiScale);
                end
                t(l) = xiScale;
            end
        end
        
        t_rm(i,:)=t;
        
        C_Cell{i}=C_rm;
        %C_Inst{i}=C_rmInst;
        
        t_teq(i,:)=t_rm(i,:)/(tEq(i));
                
    end
   
    for l=1:2 %Load Rad Mat Scaling ETCs
    %Get tOffset from statsFile if necessary
        if strcmp(tOffset,'*')==1
            disp(['** Getting tOffSet from statsFile Lag. fit for ' ...
            'tRef -> 2*tRef **'])
            tOffSet = statsFile(-1,idCell{l},runs,tref(l)*[1 2],0.5,1024,pathCell{l+1}); % dx, N kludge
        end

        %Check to see if we want to rescale
        noscaling = 0;
        if strcmp(tOffset,'noscaling')==1
            disp('** Plotting against k and removing t factor from correlator')
            tOffSet_temp=0;
            tOffSet = ones(1,nRuns)*tOffSet_temp;
            noscaling = 1;
        end
        
        %Get xiLag from statsFile 
        xiscaling = 0;
        if strcmp(tOffset,'xiscaling')==1
            disp(['** Scaling with xiLag'])
            [xiLag tStat] = statGet('xiLag',idCell{l},runs,pathCell{l+1});
            if runs(end) > 1
                xiLagAv = mean(xiLag,1);
            else
                xiLagAv = xiLag;
            end
                tOffSet_temp=0;
                tOffSet = ones(1,runs(end))*tOffSet_temp;
                xiscaling = 1;
        end

        if l==1
            [k_m,t,C,sd]=ETCload(pathCell{2},Cname,idCell{1},runs,tref(1),tOffSet,tLimit(1),tLimit(2));
        elseif l==2
            [k_r,t,C,sd]=ETCload(pathCell{3},Cname,idCell{2},runs,tref(2),tOffSet,tLimit(1),tLimit(2));
        end
        
        if (xiscaling == 1)
            %This must be done since usually tStat(end)<t(end), therefore it is
            %imposible to perform the interpolation to the lastest times.
            which=find(t<tStat(end));
            t=t(which);
            C=C(which,:);
            for q=1:size(C,1)
                xiScale = interp1(tStat,xiLagAv,t(q));
                if strcmp(Cname,'vector')~=1
                    C(q,:) = xiScale*C(q,:)/t(q);
                else
                    C(q,:) = C(q,:)*(t(q)/xiScale);
                end
                t(q) = xiScale;
            end
        end
        
        if l==1
            C_m=C;
            t_m=t;
        else
            C_r=C;
            t_r=t;
        end
        
    end

k_new=0.1:0.005:0.6;

for i=1:size(tEq,2)
    
    C_m_int=[];
    C_r_int=[];
    C_rm=[];
    C_rmInst=[];
    C_rm_int=[];
    t_new=[];
    f_t=[];
    Mean_f_t=[];
    StD_f_t=[];
    StD_f_t_temp=0;
    
    min_t=max(t_m(1),t_r(1));
    max_t=min(t_m(end),t_r(end));
    
    t_temp=t_rm(i,:);
    
    which_t=find(t_temp>min_t);
    t_new=t_rm(i,which_t);
    which_t=find(t_temp<max_t);
    t_new=t_rm(i,which_t);
    
    
     C_rm=C_Cell{i};
%     C_rmInst=C_Inst{i};
     
     C_rm=C_rm(which_t,:);
     
    %Interpolate Functions to selected k's
    for l=1:size(C_rm,1)
    
        C_rm_int(l,:)=interp1(k_rm,C_rm(l,:),k_new);
   
    end
    for l=1:size(C_m,1)
    
        C_m_temp(l,:)=interp1(k_m,C_m(l,:),k_new);
   
    end
    for l=1:size(C_r,1)
    
        C_r_temp(l,:)=interp1(k_r,C_r(l,:),k_new);
   
    end
%     for l=1:size(C_rmInst,1)
%     
%        C_rmInst_int(l,:)=interp1(kInst,C_rmInst(l,:),k_new);
%    
%     end
    %Interpolate Rad/Mat scaling functions to selected t's
    for l=1:size(C_m_temp,2)
    
        C_m_int(:,l)=interp1(t_m,C_m_temp(:,l),t_new);
   
    end
    for l=1:size(C_r_temp,2)
    
        C_r_int(:,l)=interp1(t_r,C_r_temp(:,l),t_new);
   
    end
    
    t_tEqualty{i}=t_new(:)/tEq(i);
    
    for l=1:size(C_rm_int,1)
        for j=1:size(C_rm_int,2)
            f_t(l,j)= ( C_rm_int(l,j) - C_m_int(l,j) ) / (C_r_int(l,j) - C_m_int(l,j));
        end
        
        f_joanes(i,l)=((1+(1/4)*t_new(l)/tEq(i))^-1);
        f_dani(i,l)=((1+(1/4)*t_new(l)/tEq(i))^-2);
        a(i,l)=(((2^(1/2) - 1) * (t_new(l)/tEq(i)) + 1)^2 - 1);
        f_neil(i,l)=(1+a(i,l))^(-2);
        Mean_f_t(l)=mean(f_t(l,:));
    end
    
    Mean{i}=Mean_f_t;
    
    StD_f_t=zeros(size(C_rm_int,1));
    
    for l=1:size(C_rm_int,1)
        for j=1:size(C_rm_int,2)
            StD_f_t(l) = StD_f_t(l) + (Mean_f_t(l) - f_t(l,j))^2;
        end
    end
    
    StD_f_t=StD_f_t/size(C_rm_int,2);
    StD_f_t=(StD_f_t).^(1/2);
    
    StD{i}=StD_f_t;
    
%     for j=1:size(C_rmInst_int,2)
%             f_t_Inst(j)= ( C_rmInst_int(1,j) - C_m_int(1,j) ) / (C_r_int(1,j) - C_m_int(1,j));
%     end
%     
%     Mean_f_t_Inst(i)=mean(f_t_Inst);    
%     StD_f_t_temp=0;
%     StD_f_t_Inst=0;
% 
%         for j=1:size(C_rmInst_int,2)
%             StD_f_t_temp = StD_f_t_temp + (Mean_f_t_Inst(i) - f_t_Inst(j))^2;
%         end
%     
%     StD_f_t_temp=StD_f_t_temp/size(C_rmInst_int,2);
%     StD_f_t_Inst=(StD_f_t_temp).^(1/2)
%     plotErrorBars(t_teqInst(i),Mean_f_t_Inst(i),StD_f_t_Inst,'k');hold on;

    %set(gca,'Xscale','log')
    %set(gca,'XLim',[1 10])
    
end

%figure()
for i=1:size(tEq,2)
    t_tEq=[];
    Mean_f_t=[];
    StD_f_t=[];
    t_tEq=t_tEqualty{i};
    Mean_f_t=Mean{i};
    StD_f_t=StD{i};
    plot(t_tEq,Mean_f_t); hold on;
    plotErrorBars(t_tEq(10:10:end),Mean_f_t(10:10:end),StD_f_t(10:10:end),'b');hold on;
%     %plot(t_tEq,f_dani(i,:),'k');hold on;
    %plot(t_tEq,f_joanes(i,:),'r');hold on;
    %plot(t_tEq,f_neil(i,:),'g');hold on;
    %scatter(t_teqInst(i),Mean_f_t_Inst(i),'k','filled');hold on;
    plotErrorBars(t_teqInst(i),Mean_f_t_Inst(i),StD_f_t_Inst(i),'k');
    set(gca,'Xscale','log')
end

 if(Fit==1)
     
     time=[];
     funtzion=[];
     
     %%Egin bektore bat mean eta bste bat denborarekin, dena fiteatzeko
     %%funtzioa lortzeko
     for i=1:size(tEq,2)
         tequ=t_tEqualty{i};
         mean_fun=Mean{i};
         for j=1:size(tequ,1)
            time(end+1)=tequ(j);
            funtzion(end+1)=mean_fun(j);
         end
     end

    s = fitoptions('Method','NonlinearLeastSquares','Lower',[0.1,-4],'Upper',[0.6,0]);
    f = fittype('(1 + a* x)^b','options',s);
    [c2,gof2] = fit(time',funtzion',f)
    %[c2,gof2] = fit(t_teqInst',Mean_f_t_Inst',f)
    plot(c2,'m');hold on;
 end
 
 time_teq=0.01:0.1:100;
 func_jo=(1+(1/4)*time_teq).^(-1);
 plot(time_teq,func_jo,'r');
end
   