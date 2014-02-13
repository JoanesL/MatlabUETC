%JL2013
%Usage: UETCoutputCMBeasyInterp(qmin,q,id,inPath,outPath,Ni,Evalmin,Evaldim,order,plotDotProduct)
%qmin = lowest time interval
%q  = Highest time interval
%id = 'UETCeigenScalarBeforeMat512_' id '.dat'
%           e.g: id=1
%inpath = folder which contains data
%outpath = output path relative to inpath
%Ni = size of the matrix
%Evalmin = min to plot
%Eval = number of eigenvalues to insert in CMBeasy 
%order =    -1 Reorder eigenvectors
%           -0 do not reorder
%plotDotProduct = 1 plot
%example: 11 time intervals
%
%UETCoutputCMBeasyInterp(1,11,1,pathMerged_Baseline_Tests,'InterpolatedNoOrder/',2048,1,200,1,0)

function UETCoutputCMBeasyInterp(qmin,q,id,inPath,outPath,Ni,Evalmin,Evaldim,order,plotDotProduct,rotation)

if nargin==0; 
  help UETCoutputCMBeasyInterp
  return
end

if prod(size(inPath))>0; 
  path=inPath; 
else,inPath,outPath
  if prod(size(gpath))>0
    path=gpath;
  end
end

outPath=[inPath outPath];

b=id;

for j=qmin:q
    
    disp(['Loading time interval: ' num2str(j) ])
    
    S=load([path 'UETCeigenScalarBefore_' num2str(b,'%2.2d') '_' num2str(j,'%2.2d') '.dat']);
    V=load([path 'UETCeigenVectorBefore_' num2str(b,'%2.2d') '_' num2str(j,'%2.2d') '.dat']);
    T=load([path 'UETCeigenTensorBefore_' num2str(b,'%2.2d') '_' num2str(j,'%2.2d') '.dat']);
    
    %Int=load([path 'timeIntervals_' num2str(b,'%2.2d') '.dat']);
    
    %t_frac=Int;
    
    for i=Evalmin:Evaldim
        w=i-Evalmin+1;
        s.Evector(i,:) = S(i+1,2:end);
        s.Evalue(i) = S(i+1,1);
    
        v.Evector(w,:) = V(i+1,2:end);
        v.Evalue(w) = V(i+1,1);
    
        t.Evector(i,:) = T(i+1,2:end);
        t.Evalue(i) = T(i+1,1);
    end

   if j==qmin
    ktphi=S(1,2:Ni+1);
    
    ktvt=V(1,2:end);
   end
   
   eVec_S{j}=s.Evector;
   eVec_V{j}=v.Evector;   
   eVec_T{j}=t.Evector;
      
   eVal_S{j}=s.Evalue;
   eVal_V{j}=v.Evalue;
   eVal_T{j}=t.Evalue;
    
end
   eVec{1}=eVec_S;
   eVec{2}=eVec_V;
   eVec{3}=eVec_T;
   
   eVal{1}=eVal_S;
   eVal{2}=eVal_V;
   eVal{3}=eVal_T;

    %figure()
    %a=multiPlot([1 1]);
    
    dkts=diff(ktphi(1:end));
    dktvt=diff(ktvt(1:end));        
    
    v=1:Evaldim-Evalmin+1;
    t=1:Evaldim-Evalmin+1;
    [V,T]=meshgrid(v,t);
    
    figure()
    ax=multiPlot([1 1]);
    axes(ax(1))
    
    
    %PRUEBA
    for l=q:-1:qmin
        eVal_Type=eVal{2};
        Evalplot_orig(l,:)=eVal_Type{l};
    end
    
    r=0;
    
    colour=[0 0 0]
    
    for type=1:3 %Type= S, V, T
    %for type=2
     for l=q:-1:qmin+1 %from last to first, last timestep eVectors more similar
         
         r=r+1;
         
         disp('//////////////////////')
         disp(['Interval ' num2str(l) ' vs ' num2str(l-1) ])
         disp('//////////////////////')
         
             eVec_Type=[];
             eVal_Type=[];
             Type_Evector=[];
             Type_Evalue=[];
             Type_Evector_next=[];
             Type_Evalue_next=[];
             
             
             
             eVec_Type=eVec{type};
             eVal_Type=eVal{type};

            if l==q
                Type_Evector=eVec_Type{l};
                Type_Evalue=eVal_Type{l};
                Type_Evector_next=eVec_Type{l-1};
                Type_Evalue_next=eVal_Type{l-1};
                Evalplot(l,:)=Type_Evalue(:);
            else
                Type_Evector=Type_Evector_next_Mod;
                Type_Evalue=Type_Evalue_next_Mod;
                Type_Evector_next=eVec_Type{l-1};
                Type_Evalue_next=eVal_Type{l-1};
                Evalplot(l,:)=Type_Evalue_next_Mod(:);
            end
        
            Type_Evalue_next_Mod=[];
            
            for i=1:(Evaldim+1-Evalmin) %x axis
                %k=Evalmin+i-1;
                k=i;
                for j=1:(Evaldim+1-Evalmin) % y axis
                        %s=Evalmin+j-1;
                        s=j;
                        if type==1 %Scalar
                            crossPhi(i,j)=sum(dkts.*Type_Evector(k,1:Ni-1).*(Type_Evector_next(s,1:Ni-1)));
                            crossPsi(i,j)=sum(dkts.*Type_Evector(k,Ni+1:2*Ni-1).*(Type_Evector_next(s,Ni+1:2*Ni-1)));
                            cross(i,j)=crossPhi(i,j) + crossPsi(i,j);
                        else
                            %cross(i,j)=sum(dktvt.*Type_Evector(k,1:end-1).*(Type_Evector_next(s,1:end-1)));
                            cross(i,j)=Type_Evector(k,:)*(Type_Evector_next(s,:))';
                        end
                end
            end
        
        Type_Evector_next_Mod=Type_Evector_next;
        Type_Evalue_next_Mod=Type_Evalue_next;
        
        %Reorder eVector looking at dot products
        if(order==1)
        %Always change next eVector
        up_bilatuak_index=[];
        down_bilatuak_index=[];
        for i=1:(Evaldim+1-Evalmin)
            diag=abs(cross(i,i));
            if diag<0.9
                up=0.0;
                down=0.0;
                up_index=[];
                down_index=[];
                first_up=1;
                first_down=1;
                for j=1:15
                    if ((i+j)<Evaldim-Evalmin+2 & (i-j)>0)
                    if abs(cross(i,i+j))>0.75
                        %Check and take the largest in the column
                        disp([num2str(i) ' Up found ' num2str(i+j) ' ,value ' num2str(cross(i,i+j))])
                        if first_up==1
                            first_up=0;
                            up=abs(cross(i,i+j));
                            up_index=i+j;
                        else
                            if abs(cross(i,i+j))>up;
                                up=abs(cross(i,i+j));
                                up_index=i+j;
                            end
                        end
                    elseif abs(cross(i,i-j))>0.75
                        %crossplot(i,i)=crossv(i,i-j);
                        disp([num2str(i) ' Down found ' num2str(i-j) ' ,value ' num2str(cross(i,i-j))])
                        if first_down==1
                            first_down=0;
                            down=abs(cross(i,i-j));
                            down_index=i-j;
                        else
                            if abs(cross(i,i-j))>down;
                                down=abs(cross(i,i-j));
                                down_index=i-j;
                            end
                        end
                    end
                    end
                end
                if ( up > down & up > diag)
                    disp([ num2str(i) ' Up Inserted ' num2str(up_index) ' value ' num2str(cross(i,up_index))])
                    %up_bilatuak_index(end+1)=up_index;
                    %for w=1:size(up_bilatuak_index)-1
                    %    if (up_index != up_bilatuak_index(w))
                            Type_Evector_next_Mod(i,:)=Type_Evector_next(up_index,:);
                    %        %AURRETIK=Type_Evalue_next_Mod(i)
                            Type_Evalue_next_Mod(i)=Type_Evalue_next(up_index);
                    %        %ONDOREN=Type_Evalue_next_Mod(i)
                    %    end
                    %end
                elseif ( down > up & down > diag )
                    disp([ num2str(i) ' Down Inserted ' num2str(down_index) ' value ' num2str(cross(i,down_index))])
                    %down_bilatuak_index(end+1)=up_index;
                    %for w=1:size(up_bilatuak_index)-1
                    %    if (up_index != up_bilatuak_index(w))
                            Type_Evector_next_Mod(i,:)=Type_Evector_next(down_index,:);
                    %        %AURRETIK=Type_Evalue_next_Mod(i)
                            Type_Evalue_next_Mod(i)=Type_Evalue_next(down_index);
                    %        %ONDOREN=Type_Evalue_next_Mod(i)
                    %    end
                    %end
                end
            end
        end
        if l==qmin+1
                Evalplot(l-1,:)=Type_Evalue_next_Mod(:);
        end
        end
        
%         if l==q
%             e_Vectors(l,:,:)=Type_Evector(:,:);
%         end
%             e_Vectors(l-1,:,:)=Type_Evector_next(:,:);
        %IEPE=l
        %EEEOOO=size(e_Vectors)
        
        if (rotation==1)
            e_Vectors(l,:,:)=Type_Evector_next(:,:);
            %e_vector_rot25(l,:)=Type_Evector_next(25,:);
            %e_vector_rot26(l,:)=Type_Evector_next(26,:);
            %e_vector_rot11(l,:)=Type_Evector_next(1,:);
            %e_vector_rot21(l,:)=Type_Evector_next(2,:);
        end
        
        for i=1:(Evaldim+1-Evalmin) %x axis
            %k=Evalmin+i-1;
            k=i;
            for j=1:(Evaldim+1-Evalmin) % y axis
                    %s=Evalmin+j-1;
                    s=j;
                    if type==1 %Scalar
                        crossPhiplot(i,j)=sum(dkts.*Type_Evector(k,1:Ni-1).*(Type_Evector_next_Mod(s,1:Ni-1)));
                        crossPsiplot(i,j)=sum(dkts.*Type_Evector(k,Ni+1:2*Ni-1).*(Type_Evector_next_Mod(s,Ni+1:2*Ni-1)));
                        crossplot(i,j)=crossPhiplot(i,j) + crossPsiplot(i,j);
                    else
                        %crossplot(i,j)=sum(dktvt.*Type_Evector(k,1:end-1).*(Type_Evector_next_Mod(s,1:end-1)));
                        crossplot(i,j)=Type_Evector(k,:)*(Type_Evector_next_Mod(s,:)');
                    end
            end
            if crossplot(i,i)<0 %Relative sign
                 %disp([ num2str(type) ' sign ' num2str(i)])
                 Type_Evector_next_Mod(k,:)=-Type_Evector_next_Mod(k,:);
                 crossplot(i,i)=-crossplot(i,i);
            end
        end
        %crosss=crosssphi+crossspsi;        
        if type==1
            Cname='Scalar';
        elseif type==2
            Cname='Vector';
        elseif type==3
            Cname='Tensor';
        end
                
%         %plot(ktphi, Type_Evector(36,1:Ni), 'k',ktpsi, Type_Evector(36,Ni+1:2*Ni), 'b'); hold on;
%         colour(2)=colour(2)+0.03;
%         axes(ax(1))
%         plot(ktvt, Type_Evector_next(20,1:Ni),'Color',colour);hold on;
%         %set(gca,'Color',colour);hold on;
%         titleplot=[Cname ' e-vector evolution ' num2str(25) ];
%         title(titleplot ,'FontWeight','bold')
%         set(gca,'Xscale','log')
        
        if plotDotProduct==1
        figure()
        plotc(V,T,abs(cross),1,Evaldim-Evalmin+1,100,0) 
        titleplot=['Cross ' Cname ' ' num2str(l) ' vs ' num2str(l-1) ];
        title(titleplot ,'FontWeight','bold')
        
        if order==1
        figure()
        plotc(V,T,crossplot,1,Evaldim-Evalmin+1,100,0) 
        titleplot=['Cross ' Cname ' Modified ' num2str(l) ' vs ' num2str(l-1) ];
        title(titleplot ,'FontWeight','bold')
        end
        end
        
        if l==q
            if type==1 %S
                eigenoutput(Cname,ktphi,Type_Evalue_next,Type_Evector_next,Ni,Evaldim,1,outPath,b,l);
            else
                eigenoutput(Cname,ktvt,Type_Evalue_next,Type_Evector_next,Ni,Evaldim,0,outPath,b,l);
            end
        end
        if type==1 %S
            eigenoutput(Cname,ktphi,Type_Evalue_next_Mod,Type_Evector_next_Mod,Ni,Evaldim,1,outPath,b,l-1);
        else
            eigenoutput(Cname,ktvt,Type_Evalue_next_Mod,Type_Evector_next_Mod,Ni,Evaldim,0,outPath,b,l-1);
        end
        
     end
     
     %e_Vector50_temp_esk(1:Ni)=e_Vectors(1,50,:);
     %e_Vector51_temp_esk(1:Ni)=e_Vectors(1,51,:);
     
     
%      for i=1:q
%      %e_Vector50_temp_esk(1:Ni)=e_Vectors(i,50,:);
%      %e_Vector51_temp_esk(1:Ni)=e_Vectors(i,51,:);
%         for j=1:Evaldim
%             e_Vector_temp(1:Ni)=e_Vectors(i,j,:);
%             dot_pro50(j,i)=e_Vector50_temp_esk*e_Vector_temp';
%             dot_pro51(j,i)=e_Vector51_temp_esk*e_Vector_temp';
%         end
%      end
%      
%      figure()
%      plot(abs(dot_pro50))
%      grid on;
%      set(gca,'XLim',[40 60])
%      set(gca,'YLim',[0 1.05])
%      xlabel('e-Vector index')
%      ylabel('Dot Product')
%      title('LEFT e-Vector 50 Dot Product','FontWeight','bold') 
%      
%      figure()
%      plot(abs(dot_pro51))
%      grid on;
%      set(gca,'XLim',[40 60])
%      set(gca,'YLim',[0 1.05])
%      xlabel('e-Vector index')
%      ylabel('Dot Product')
%      title('LEFT e-Vector 51 Dot Product','FontWeight','bold')
     
     
     if (rotation==1)
         
         t_i=[8.967 11.27 14.49];
         for i=1:size(t_i)-1
             t_bar(i)=(t_i(i)+t_i(i+1))/2;
         end
         
         t=linspace(t_i(1), t_i(3), 25);
         color_frac=linspace(0,1,25);
         
         JUASS=size(e_Vectors)
                  
         theta_t=pi * (t-t_i(1))/(t_i(3)-t_i(1));
              
         for j=1:50 %Number of e-vectors.
            for i=1:size(t,2)
                e_Vectors_Rot(j,i,:)=cos(theta_t(i)/2)*e_Vectors(9,j,:)+sin(theta_t(i)/2)*e_Vectors(10,j,:);
            %e_Vector25(i,:)=cos(theta_t(i)/2)*e_vector_rot25(9,:)+sin(theta_t(i)/2)*e_vector_rot25(10,:);
            %e_Vector26(i,:)=cos(theta_t(i)/2)*e_vector_rot26(9,:)+sin(theta_t(i)/2)*e_vector_rot26(10,:);
            %e_Vector1(i,:)=cos(theta_t(i)/2)*e_vector_rot11(9,:)+sin(theta_t(i)/2)*e_vector_rot11(10,:);
            %e_Vector2(i,:)=cos(theta_t(i)/2)*e_vector_rot21(9,:)+sin(theta_t(i)/2)*e_vector_rot21(10,:);
            %dot_pro(i)=e_Vector25(i,:)*e_Vector2(i,:)';
            %dot_pro2(i)=e_Vector1(i,:)*e_Vector25(i,:)';
            end
         end
         
         for j=1:50
         for i=1:size(t,2)
             e_Vector25_temp(1:Ni)=e_Vectors_Rot(10,i,:);
             e_Vectors_temp(1:Ni)=e_Vectors_Rot(j,i,:);
             EEOOO=size(e_Vectors_temp)
             dot_pro(j,i)=e_Vector25_temp*e_Vectors_temp';
         end
         end
         
         EEEOOO=size(e_Vectors_Rot)
         sizeDOTPRO=size(dot_pro)
         
         colour=zeros(25,3);
         
         for i=1:25
             colour(i,:)=[0 0 0+color_frac(i)];
         end
         
         
        figure() 
        %plot(ktvt,e_vector_rot1(9,:),'b',ktvt,e_vector_rot1(10,:),'k');
        %set(gca,'Xscale','log')
        plot(dot_pro);
        %set(gca,'Xscale','log')
        %set(gca, 'Xlim', [t_i(1) t_i(3)])
        
        %figure()
        %set(0,'DefaultAxesColorOrder',colour)
        %for i=1:size(t,2)
        %plot(ktvt,e_Vector26);hold on;
        %end
        
        %set(gca,'Color',colour);hold on;
        titleplot=[Cname ' e-vector rotation ' num2str(25) ];
        title(titleplot ,'FontWeight','bold')
        set(gca,'Xscale','log')
         
     end
     
     
     %if type==2
     %    q_index=q:-1:qmin+1;
         %q_index=qmin:1:q-1
     %    figure()
     %    plot(q_index,abs(eVal22),'b',q_index,abs(eVal23),'k');hold on;
     %end
    end
    
    figure()
    plot(abs(Evalplot_orig(:,15:25)),'linewidth',1.2);hold on; grid on;
        set(gca,'XLim',[qmin q])
        set(gca,'XTick',[qmin:1:q])
        xlabel('Interval index')
        ylabel('Evalue')
        title('Evalue evolution Before reordering','FontWeight','bold')
    
    
    figure()
        plot(abs(Evalplot(:,15:25)),'linewidth',1.2);hold on; grid on;
        set(gca,'XLim',[qmin q])
        set(gca,'XTick',[qmin:1:q])
        xlabel('Interval index')
        ylabel('Evalue')
        title('Evalue evolution After reordering','FontWeight','bold')    
    
    
    %figure()
    %for n=15:30
    %    plot(qmin:1:q,abs(Evalplot(:,n)));hold on
    %end
    %figure()
    %for n=15:30
    %    plot(qmin:1:q,Evalplot(:,n));hold on
    %end
    
	    %disp('Lets check if the order was correctly set'  


end

function eigenoutput(Cname,kt,eValue,eVector,Ni,Evaldim,scalar,outpath,id,q)

if scalar==1
    ktI=logspace(-3,log10(kt(end)),8*Ni);
else
    ktI=logspace(-3,log10(kt(end)),4*Ni);
end

if scalar==1
    for i=1:Evaldim 
        eVectorPhiI(i,:)=interp1([0 kt],[0 eVector(i,1:Ni)], ktI, 'pchip');
        eVectorPsiI(i,:)=interp1([0 kt],[0 eVector(i,Ni+1:end)], ktI, 'pchip');
    end
else
    for i=1:Evaldim 
        eVectorI(i,:)=interp1([0 kt],[0 eVector(i,:)], ktI, 'pchip');
    end
end


if scalar==1
    M=[2*max(size(ktI)) ktI ktI; (abs(eValue))' eVectorPhiI eVectorPsiI];
else
    M=[max(size(ktI)) ktI; (abs(eValue))' eVectorI];
end

save([outpath 'UETCeigen' Cname '_' num2str(id,'%2.2d') '_' num2str(q,'%2.2d') '.dat'],'M','-ascii')

disp('UETC interpolated eigenvectors saved...')
end
