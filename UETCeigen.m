%UETC eigenvector decomposition function for UETC.hpp data 2006-2009
%
%Change history: April 2009 - added highKtMode input and extrapolation at high kt 
%
%Usage: UETCeigen(id,runs,tRef,tOffSet,inPath,outPath,era,number,bootstrap,Ni,before,rMax,plateauMode,highKtMode,plotLevel)
%
%   id = ID string between 'UETCscalar11_' and before '.dat'
%          eg. to load statsFile_6L01 ID is '6L%2' with run=1
%    run = realizations(s) to include
%   tRef = UETC reference time
%tOffset = time when xi=0, if '*' get from statsFile (Lag & tRef<t<2*tRef)
%          (independent offsets for each realization is used, 
%           but a common one (ie. single number) can be specified too)
%
%Optional parameters are:
%
%   inPath = path to input UETCs files, including final '/'
%            (if ommited or '' gets path from gpath global variable)
%  outPath = path (relative to inPath) for outputting UETCeigen files
%      era = either 'Mat' or 'Rad' (simply for output file names)
%   number = ID number in output filenames, default 0
%bootstrap = 0 means no, use average 
%            1 means use average of random 
%            realizations combinations (without removal)
%       Ni = (interpolated) matrix size (default 512)
%		Before= 1 output eigenvectors before interpolation
%     rMax = zero data above rMax
%plateauMode = 0 means use lowest kt point (point 1)
%              1 means use mean of points 2 and 3 (and replace 1)
%              (default is 1)
%highKtMode = 0 means use raw data for high Kt
%             1 means extrapolate/correct as power-law at high Kt
%plotLevel = 0 for zero plots
%            1 for more
%            2 for even more

function UETCeigen(id,runs,tRef,tOffSet,inPath,outPath,era,number,bootstrap,Ni,before,rMax,plateauMode,highKtMode,plotLevel)

if nargin==0; 
  help UETCeigen
  return
end

global gpath

if ~exist('inPath','var'); inPath=''; end 

if prod(size(inPath))>0; 
  path=inPath; 
else,inPath,outPath,era
  if prod(size(gpath))>0
    path=gpath;
  else
    disp(['Please set gpath global variable to default path'...
	  ' or specify path in fucntion call'])
    return
  end
end

if ~exist('outPath','var'); outPath=''; end
if ~exist('bootstrap','var'); bootstrap=0; end
if ~exist('Ni','var'); Ni=512; end
if ~exist('rMax','var'); rMax=99; end
if rMax==0; rMax=99; end
if ~exist('plateauMode','var'); plateauMode=1; end
if ~exist('plotLevel','var'); plotLevel=1; end
if ~exist('highKtMode','var'); highKtMode=1; end
  
%Perform bootstrapping over realization numbers if bootstrap set
if bootstrap==1
   runs=runs( ceil( size(runs,2)*rand(size(runs)) ) )
elseif bootstrap ~= 0
  disp('Bootstrap input should be 0 or 1!')
  return
end

%Get tOffset from statsFile if necessary
if strcmp(tOffSet,'*')==1
  disp(['** Getting tOffSet from statsFile Lag. fit for ' ...
	'tRef -> 2*tRef **'])
  tOffSet = statsFile(-1,id,runs,tRef*[1 2]);
  disp(['** Result: tOffSet = ' num2str(tOffSet) ' **'])
end

%Expand tOffSet to vector if only one given for multiple runs
if prod(size(tOffSet))==1; tOffSet=tOffSet*ones(size(runs)); end

%Get interpolation coordinates (x=kt)
disp(['Matrix size: ' num2str(Ni) 'x' num2str(Ni)])

%x=linspace(0.4,200,Ni);
%disp(['Linear up to 200 (2005 default)'])

%x=linspace(0.4,1500,Ni);
%disp(['Linear up to 1500 (2009 default)'])

%x=logspace(log10(0.1),log10(1500),Ni);
%disp(['Log up to 1500'])


%Joanes Lizarraga (28/11/13)
%Hardcoded for Merged UETCs
%x=linspace(0.7,2350,Ni);
%x=linspace(0.1,2350,Ni);
x=logspace(log10(0.1),log10(2350),Ni);
%s=0
%x=linspace(0.75,2000,Ni);

%Log for low kt, linear for high kt
%x1=logspace(log10(0.1),log10(3),100);
%x2=linspace(3.05,2350,Ni-100);

%x=[x1 x2];



%x1=linspace(0.4,200,65);
%x2=logspace(log10(x1(end)),log10(2000),Ni-64);
%x=[x1 x2(2:end)];
%disp(['Linear up to 200 with Ni-64 points, then logspaced to 1500'])

%x=logspace(log10(0.4),log10(200),Ni);
%x=logspace(log10(0.4),log10(8),Ni/2+1); x=[x(1:end-1) linspace(8,200,Ni/2)];

disp(['Matrix kt from:' num2str(min(x)) ' to ' num2str(max(x)) ' (hardcoded)'])

[X,Xdash]=meshgrid(x,x);

for i=1:prod(size(runs))
  
   %LOAD CORRELATORS
   [kt,r,s.C11]=UETCload(path,'scalar11',id,runs(i),tRef,tOffSet(i),0); 
   [kt,r,s.C12]=UETCload(path,'scalar12',id,runs(i),tRef,tOffSet(i),0);
   [kt,r,s.C21]=UETCload(path,'scalar21',id,runs(i),tRef,tOffSet(i),0);
   [kt,r,s.C22]=UETCload(path,'scalar22',id,runs(i),tRef,tOffSet(i),0);
   [kt,r,v.C]=UETCload(path,'vector',id,runs(i),tRef,tOffSet(i),0);
   [kt,r,t.C]=UETCload(path,'tensor',id,runs(i),tRef,tOffSet(i),0);
      
   ktMin=kt(1)
   ktMax=kt(end)
   
   %DON'T CORRECT TENSOR NORMALIZATION ERROR
   %t.C=t.C/2;
   disp('** Note that this version assumes correct tensor normalization **')
   
   %WORRY ABOUT ktmax and max(x)
   disp(['Max x: ' num2str(x(end)) ' cf. max kt: ' num2str(kt(end))])

   %ZERO DATA AT LARGE RATIOS (if needed)
   disp(['Max r in data: ' num2str(max(r))])
   if max(r)>rMax
      disp(['Cutting max(r) from: ' num2str(max(r)) 'to: ' num2str(max(r(r<=rMax)))])
      s.C11(r>rMax,:)=0;
      s.C12(r>rMax,:)=0;
      s.C21(r>rMax,:)=0;
      s.C22(r>rMax,:)=0;
      v.C(r>rMax,:)=0;
      t.C(r>rMax,:)=0;
   end
   
   %EXTRAPOLATE/CORRECT DATA AT HIGH KT
   %(should be done first since uses ETCs data)
   
   if highKtMode==1
     [temp,temp,ETC1]=ETCload(path,'scalar11',id,runs(i),tRef,tOffSet(i),0,999);
     s.C11=UETCextrap(s.C11,kt,r,ETC1);

     [temp,temp,ETC2]=ETCload(path,'scalar22',id,runs(i),tRef,tOffSet(i),0,999);
     s.C22=UETCextrap(s.C22,kt,r,ETC2);

     s.C12=UETCextrap(s.C12,kt,r,ETC1,ETC2);
     s.C21=UETCextrap(s.C21,kt,r,ETC2,ETC1);

     [temp,temp,ETC1]=ETCload(path,'vector',id,runs(i),tRef,tOffSet(i),0,999);
     v.C=UETCextrap(v.C,kt,r,ETC1);

     [temp,temp,ETC1]=ETCload(path,'tensor',id,runs(i),tRef,tOffSet(i),0,999);
     t.C=UETCextrap(t.C,kt,r,ETC1);
   end

   %ADD EXTRA DATA AT LOW kt
   if plateauMode==0
     disp('Plateau mode: as lowest kt point')
     ktextra=min(x):min(x):kt(1);
     kt=[ktextra kt(1:end)];
     s.C11=[repmat(s.C11(:,1),size(ktextra)) s.C11(:,1:end)];
     s.C12=[repmat(s.C12(:,1),size(ktextra)) s.C12(:,1:end)];
     s.C21=[repmat(s.C21(:,1),size(ktextra)) s.C21(:,1:end)];
     s.C22=[repmat(s.C22(:,1),size(ktextra)) s.C22(:,1:end)];
     v.C=[repmat(v.C(:,1),size(ktextra)) v.C(:,1:end)];
     t.C=[repmat(t.C(:,1),size(ktextra)) t.C(:,1:end)];
   elseif plateauMode==1
     disp('Plateau mode: as mean of low kt points 2 and 3')
     ktextra=min(x):min(x):kt(3);
     kt=[ktextra kt(4:end)];
     s.C11=[repmat(mean(s.C11(:,2:3),2),size(ktextra)) s.C11(:,4:end)];
     s.C12=[repmat(mean(s.C12(:,2:3),2),size(ktextra)) s.C12(:,4:end)];
     s.C21=[repmat(mean(s.C21(:,2:3),2),size(ktextra)) s.C21(:,4:end)];
     s.C22=[repmat(mean(s.C22(:,2:3),2),size(ktextra)) s.C22(:,4:end)];
     v.C=[repmat(mean(v.C(:,2:3),2),size(ktextra)) v.C(:,4:end)];
     t.C=[repmat(mean(t.C(:,2:3),2),size(ktextra)) t.C(:,4:end)];
   else
   end

   %GENERATE RAW COORDINATE ARRAYS
   Nk=size(kt,2);
   Nr=size(r,2);
   [KT,R]=meshgrid(kt,r);
   
   %MAP TO KT-KT' PLANE
   [KTmap,Rmap]=UETCktr2xx(KT,R,X,Xdash);
   if i==1;        
      s.M=UETCmatrixScalar(KTmap,Rmap,s.C11,s.C12,s.C21,s.C22);
      v.M=UETCmatrix(KTmap,Rmap,v.C);
      t.M=UETCmatrix(KTmap,Rmap,t.C);
   else
      s.M=s.M+UETCmatrixScalar(KTmap,Rmap,s.C11,s.C12,s.C21,s.C22);
      v.M=v.M+UETCmatrix(KTmap,Rmap,v.C);
      t.M=t.M+UETCmatrix(KTmap,Rmap,t.C) ;       
   end    
end

%READY FIGURE (now all files are found)
if plotLevel>3
  figure
  a=multiPlot([3 3]);
end

%AVERAGE MATRICES
disp('AVERAGING MATRICES...')
s.M=s.M/prod(size(runs));
v.M=v.M/prod(size(runs));
t.M=t.M/prod(size(runs));

if plotLevel>3
   %PLOT s.M MATRIX
   axes(a(1))
   pcolor(log10(abs(s.M)))
   xlabel('x index')
   ylabel('x'' index')
   view([90 90])
   axis tight;
   shading interp
   
   %PLOT v.M MATRIX
   axes(a(2))
   pcolor(log10(abs(v.M)))
   xlabel('x')
   ylabel('x''')
   view([90 90])
   axis tight
   shading interp
   
   %PLOT t.M MATRIX
   axes(a(3))
   pcolor(log10(abs(t.M)))
   xlabel('x')
   ylabel('x''')
   view([90 90])
   axis tight
   shading interp
end

%CALCULATE EIGENVALUES AND VECTORS
disp('CALCULATING EIGENVALUES AND EIGENVECTORS...')
[s.Evector,s.Evalue]=eigen(s.M);
[v.Evector,v.Evalue]=eigen(v.M);
[t.Evector,t.Evalue]=eigen(t.M);



if before==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%JOJOJOJO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save eigenvecotors: usual format and format for Reconstruction
fileeigensca=[path outPath 'UETCeigenScalarBefore' era '_' num2str(number,'%2.2d') '.dat'];
fileeigenvec=[path outPath 'UETCeigenVectorBefore' era '_' num2str(number,'%2.2d') '.dat'];
fileeigenten=[path outPath 'UETCeigenTensorBefore' era '_' num2str(number,'%2.2d') '.dat'];




A=[2*max(size(x)) x x; s.Evalue (s.Evector)'];
B=[max(size(x)) x; v.Evalue (v.Evector)'];
C=[max(size(x)) x; t.Evalue (t.Evector)'];


save(fileeigensca,'-ascii','A');
save(fileeigenvec,'-ascii','B');
save(fileeigenten,'-ascii','C');
disp('Saving eigenvalues and eigenvectors...')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%JOJOJOJO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %figure()
    %semilogx(x, s.Evector(1:Ni,1),'r', x, s.Evector(1:Ni,2),'g', x, s.Evector(1:Ni,3),'b',x, s.Evector(1:Ni,1)+s.Evector(1:Ni,2)+s.Evector(1:Ni,3),'k'); hold on %JO changed 3 first eigenvectors
  %semilogx(x, s.Evector(Ni+1:end,1),'--r', x, s.Evector(Ni+1:end,2),'--g', x, s.Evector(Ni+1:end,3),'--b'); hold on
  %xlabel('index')
  %ylabel('Eigenvector value')
  %axis tight

if plotLevel==1
  disp('CREATING PLOTS...')

  %PLOT SCALAR EIGENVALUES
  axes(a(4))
  plot(1:2*Ni, s.Evalue,'b', 1:2*Ni,abs(s.Evalue),'r', [1 2*Ni],[0 0],'k'); hold on
  xlabel('Eigenvalue number')
  ylabel('Eigenvalue')
  axis tight
  set(gca,'Xlim',[1 50])

  %PLOT VECTOR EIGENVALUES
  axes(a(5))
  plot(1:Ni, v.Evalue,'b', 1:Ni,abs(v.Evalue),'r', [1 Ni],[0 0],'k')
  xlabel('Eigenvalue number')
  ylabel('Eigenvalue')
  axis tight
  set(gca,'Xlim',[1 50])
  
  %PLOT TENSOR EIGENVALUES
  axes(a(6))
  plot(1:Ni, t.Evalue,'b', 1:Ni,abs(t.Evalue),'r', [1 Ni],[0 0],'k')
  xlabel('Eigenvalue number')
  ylabel('Eigenvalue')
  axis tight
  set(gca,'Xlim',[1 50])
 
 
  %PLOT FIRST 3 SCALAR EIGENVECTORS
  axes(a(7))
  semilogx(x, s.Evector(1:Ni,1),'r', x, s.Evector(1:Ni,2),'g', x, s.Evector(1:Ni,3),'b'); hold on %JO changed 3 first eigenvectors
  semilogx(x, s.Evector(Ni+1:end,1),'--r', x, s.Evector(Ni+1:end,2),'--g', x, s.Evector(Ni+1:end,3),'--b'); hold on
  xlabel('index')
  ylabel('Eigenvector value')
  axis tight

  %PLOT FIRST 3 VECTOR EIGENVECTORS
  axes(a(8))
  semilogx(x, v.Evector(:,1),'r', x, v.Evector(:,2),'g', x, v.Evector(:,3),'b'); hold on
  xlabel('index')
  ylabel('Eigenvector value')
  axis tight

  %PLOT FIRST 3 TENSOR EIGENVECTORS
  axes(a(9))
  semilogx(x, t.Evector(:,1),'r', x, t.Evector(:,2),'g', x, t.Evector(:,3),'b'); hold on
  xlabel('index')
  ylabel('Eigenvector value')
  axis tight
end
  
%WRITE OUTPUT
%if prod(size(outPath))>0
   disp('WRITING OUTPUT...')
   output([path outPath 'UETCeigenScalar' era '.dat'],x,s.Evalue,s.Evector,1,Ni,0);
   output([path outPath 'UETCeigenVector' era '.dat'],x,v.Evalue,v.Evector,0,Ni,1);
   output([path outPath 'UETCeigenTensor' era '.dat'],x,t.Evalue,t.Evector,0,Ni,0);
%end

if plotLevel==2
   
   disp('RECONSTRUCTING...')

   %RECONSTRUCT ETCs 
   figure
   aa=multiPlot([3 2]);
   delete(aa(6)); aa=aa(1:5);
   
   ETC11=zeros(Ni,1);
   ETC22=ETC11; 
   ETC12=ETC11;
   ETCv=ETC11;
   ETCt=ETC11;
   for i=1:2*Ni
      ETC11=ETC11+s.Evalue(i)*s.Evector(1:Ni,i).^2;
      ETC22=ETC22+s.Evalue(i)*s.Evector(Ni+1:end,i).^2;
      ETC12=ETC12+s.Evalue(i)*s.Evector(1:Ni,i).*s.Evector(Ni+1:end,i);
      if i==1 | i==2 | i==3 | i==100 | i==200 | i==2*Ni
         if i==2*Ni; color='k'; elseif i==1; color='b'; elseif i==2; color='r';
         elseif i==3; color='g'; elseif i==100; color='c'; else color='m'; end
         axes(aa(1)); plot(x,ETC11,color); hold on; set(gca,'XScale','log');
         axes(aa(2)); plot(x,ETC12,color); hold on; set(gca,'XScale','log');
         axes(aa(3)); plot(x,ETC22,color); hold on; set(gca,'XScale','log');
      end
   end
   for i=1:Ni
      ETCv=ETCv+v.Evalue(i)*v.Evector(:,i).^2;
      ETCt=ETCt+t.Evalue(i)*t.Evector(:,i).^2;
      if i==1 | i==2 | i==3 | i==100 | i==200 | i==Ni
         if i==Ni; color='k'; elseif i==1; color='b'; elseif i==2; color='r';
         elseif i==3; color='g'; elseif i==100; color='c'; else color='m'; end
         axes(aa(4)); plot(x,ETCv,color); hold on; set(gca,'XScale','log');
         axes(aa(5)); plot(x,ETCt,color); hold on; set(gca,'XScale','log');
      end
   end
   for i=1:5; axes(aa(i)); axis tight; end
   multiPlotZoom(aa)
  
   %RELOAD RAW UETCs (to check reconstructions)
   [kt2,r2,UETC11]=UETCload(path,'scalar11',id,runs(1),tRef,tOffSet(1),0); 
   [kt2,r2,UETC12]=UETCload(path,'scalar12',id,runs(1),tRef,tOffSet(1),0); 
   [kt2,r2,UETC21]=UETCload(path,'scalar21',id,runs(1),tRef,tOffSet(1),0); 
   [kt2,r2,UETC22]=UETCload(path,'scalar22',id,runs(1),tRef,tOffSet(1),0); 
   [kt2,r2,UETCv]=UETCload(path,'vector',id,runs(1),tRef,tOffSet(1),0); 
   [kt2,r2,UETCt]=UETCload(path,'tensor',id,runs(1),tRef,tOffSet(1),0); 
   
   %PLOT ETC RECONSTRUCTIONS
   figure
   aa=multiPlot([3 2]);
   delete(aa(6)); aa=aa(1:5);
   
   axes(aa(1))
   plot(kt2(1:end),UETC11(1,:),'o-k'); hold on
   plot(x,ETC11,'.-')
   set(gca,'XScale','log')
   xlabel('kt')
   ylabel('ETC C_1_1')
   axis tight
   
   axes(aa(2))
   plot(kt2(1:end),UETC12(1,:),'o-k'); hold on
   plot(x,ETC12,'.-')
   set(gca,'XScale','log')
   xlabel('kt')
   ylabel('ETC C_1_2')
   axis tight
   
   axes(aa(3))
   plot(kt2(1:end),UETC22(1,:),'o-k'); hold on
   plot(x,ETC22,'.-')
   set(gca,'XScale','log')
   xlabel('kt')
   ylabel('ETC C_2_2')
   axis tight
   
   axes(aa(4))
   plot(kt2(1:end),UETCv(1,:),'o-k'); hold on
   plot(x,ETCv,'.-')
   set(gca,'XScale','log')
   xlabel('kt')
   ylabel('ETC C^V')
   axis tight
   
   axes(aa(5))
   plot(kt2(1:end),UETCt(1,:),'o-k'); hold on
   plot(x,ETCt,'.-')
   set(gca,'XScale','log')
   xlabel('kt')
   ylabel('ETC C^T')
   axis tight
   
   multiPlotZoom(aa)

end

if plotLevel==3
   
   %PLOT ORIGINAL UETCs
   %figure
   %aaa=multiPlot([3 2]);
   
   %axes(aaa(1))
   %[KT,R]=meshgrid(kt2(kt2<max(x)),r2);
   %surf(KT.*sqrt(R),R,UETC11(:,kt2<max(x)))
   %hold on
   
   %axes(aaa(2))
   %[KT,R]=meshgrid(kt2(kt2<max(x)),r2);
   %surf(KT.*sqrt(R),R,UETC12(:,kt2<max(x)))
   %hold on
   
   %axes(aaa(3))
   %[KT,R]=meshgrid(kt2(kt2<max(x)),r2);
   %surf(KT.*sqrt(R),R,UETC21(:,kt2<max(x)))
   %hold on
   
   %axes(aaa(4))
   %[KT,R]=meshgrid(kt2(kt2<max(x)),r2);
   %surf(KT.*sqrt(R),R,UETC22(:,kt2<max(x)))
   %hold on
   
   %axes(aaa(5))
   %[KT,R]=meshgrid(kt2(kt2<max(x)),r2);
   %surf(KT.*sqrt(R),R,UETCv(:,kt2<max(x)))
   %hold on
   
   %axes(aaa(6))
   %[KT,R]=meshgrid(kt2(kt2<max(x)),r2);
   %surf(KT.*sqrt(R),R,UETCt(:,kt2<max(x)))
   %hold on
   
   %RESCONSTRUCT UETCs
   %UETCs=zeros(2*Ni,2*Ni);
   %for j=1:2*Ni
   %   for k=1:2*Ni
   %      UETCs(j,k)=(s.Evector(j,:).*s.Evector(k,:))*abs(s.Evalue);
   %   end
   %end
   
   %UETCs=s.M./UETCs;
   
   %UETC11=UETCs(1:Ni,1:Ni);
   %UETC12=UETCs(1:Ni,Ni+1:end);
   %UETC21=UETCs(Ni+1:end,1:Ni);
   %UETC22=UETCs(Ni+1:end,Ni+1:end);
   %clear UETCs;
   
   
   UETCv_OK=zeros(Ni,Ni);
   %UETCt=zeros(Ni,Ni);
   for j=1:Ni
      for k=1:Ni
         UETCv_OK(j,k)=(v.Evector(j,:).*v.Evector(k,:))*v.Evalue;
         %UETCt(j,k)=(t.Evector(j,:).*t.Evector(k,:))*abs(t.Evalue);
      end
   end
   
    for i=1:size(v.Evalue)
       if(v.Evalue(i)<0)
           disp(['Negative! ' num2str(i)])
           v.Evalue(i)=0.00001*v.Evalue(i);
       end
    end
    UETCv=zeros(Ni,Ni);
   %UETCt=zeros(Ni,Ni);
   for j=1:Ni
      for k=1:Ni
         UETCv(j,k)=(v.Evector(j,:).*v.Evector(k,:))*v.Evalue;
         %UETCt(j,k)=(t.Evector(j,:).*t.Evector(k,:))*abs(t.Evalue);
      end
   end
   
   ABS_V=UETCv-v.M;
   OK_V=UETCv_OK-v.M;
   
   [X,Xdash]=meshgrid(x,x);
   KT=sqrt(X.*Xdash);
   R=X./Xdash;
   
   %axes(aaa(1))
   %surf(KT,R,UETC11)%./(R<=1 & R>0.5/max(r)))
   %xlabel('k(tt'')^0^.^5')
   %ylabel('t''/t')
   %axis tight
   %set(gca,'YLim',[0.5./max(r) max(r)])
   %set(gca,'XLim',[min(x) max(x)*max(r)])
   %view([0 90])
   %view([-141 36])
   %shading interp
   %set(gca,'XScale','log')
   %set(gca,'YScale','log')
   
   %axes(aaa(2))
   %surf(KT,R,UETC12)%./(R<=1 & R>0.5/max(r)))
   %xlabel('k(tt'')^0^.^5')
   %ylabel('t''/t')
   %axis tight
   %set(gca,'YLim',[0.5./max(r) max(r)])
   %set(gca,'XLim',[min(x) max(x)*max(r)])
   %view([0 90])
   %view([-141 36])
   %shading interp
   %set(gca,'XScale','log')
   %set(gca,'YScale','log')
   
   %axes(aaa(3))
   %surf(KT,R,UETC21)%./(R<=1 & R>0.5/max(r)))
   %xlabel('k(tt'')^0^.^5')
   %ylabel('t''/t')
   %axis tight
   %set(gca,'YLim',[1/max(r) max(r)])%[0.5./max(r) max(r)])
   %set(gca,'XLim',[min(x) max(x)*max(r)])
   %view([0 90])
   %view([-141 36])
   %shading interp
   %set(gca,'XScale','log')
   %set(gca,'YScale','log')
   
   %axes(aaa(4))
   %surf(KT,R,UETC22)%./(R<=1 & R>0.5/max(r)))
   %xlabel('k(tt'')^0^.^5')
   %ylabel('t''/t')
   %axis tight
   %set(gca,'YLim',[1/max(r) max(r)])%[0.5./max(r) max(r)])
   %set(gca,'XLim',[min(x) max(x)*max(r)])
   %view([0 90])
   %view([-141 36])
   %shading interp
   %set(gca,'XScale','log')
   %set(gca,'YScale','log')
   
   %axes(aaa(5))
   figure()
   colour=colourMap(ABS_V);
   %surf(KT,R,UETCv)%./(R<=1 & R>0.5/max(r)))
   surf(KT,R,ABS_V,colour,'FaceLighting','phong','FaceColor','interp','EdgeColor','flat');
   xlabel('k(tt'')^0^.^5')
   ylabel('t''/t')
   axis tight
   set(gca,'YLim',[1/max(r) max(r)])%[0.5./max(r) max(r)])
   set(gca,'XLim',[min(x) max(x)*max(r)])
   %view([0 90])
   view([-141 36])
   shading interp
   set(gca,'XScale','log')
   set(gca,'YScale','log')
   
   figure()
   colour=colourMap(OK_V);
   %surf(KT,R,UETCv)%./(R<=1 & R>0.5/max(r)))
   surf(KT,R,OK_V,colour,'FaceLighting','phong','FaceColor','interp','EdgeColor','flat');
   xlabel('k(tt'')^0^.^5')
   ylabel('t''/t')
   axis tight
   set(gca,'YLim',[1/max(r) max(r)])%[0.5./max(r) max(r)])
   set(gca,'XLim',[min(x) max(x)*max(r)])
   %view([0 90])
   view([-141 36])
   shading interp
   set(gca,'XScale','log')
   set(gca,'YScale','log')
   
   %axes(aaa(6))
   %surf(KT,R,UETCt)%./(R<=1 & R>0.5/max(r)))
   %xlabel('k(tt'')^0^.^5')
   %ylabel('t''/t')
   %axis tight
   %set(gca,'YLim',[1/max(r) max(r)])%[0.5./max(r) max(r)])
   %set(gca,'XLim',[min(x) max(x)*max(r)])
   %view([0 90])
   %view([-141 36])
   %shading interp
   %set(gca,'XScale','log')
   %set(gca,'YScale','log')
end

runs

%ACTIVATE MULTIPLOT ZOOM
%multiPlotZoom(a)


%=============================================================
%=============================================================

function [vector,value]=eigen(C)

N=size(C,1);

%Use eig to find eigenvalues and vectors
tic
[vector,value]=eig(C);
t=toc;
disp(['Eigenvectors determined in ' num2str(t) 's'])
value=diag(value);

%Sort in descending eigenvalue (magnitude... if desired)
[dummy,index]=sort(-abs(value));
%[dummy,index]=sort(-value);
value=value(index);
vector=vector(:,index);

%Ensure first point in eigenvector (low kt) is postive
%this is a sign convention so that matter and radiation eras have same
%sign... I hope!
%vector=vector.*repmat(sign(vector(1,:)),N,1);

%Or use the mean must be postive
vector=vector.*repmat(sign(mean(vector,1)),N,1);

%=============================================================
%=============================================================

function [KTmap,Rmap]=UETCktr2xx(KT,R,X,Xdash)

Ni=size(X,1);
KTmap=zeros(Ni,Ni);
Rmap=zeros(Ni,Ni);

%Loop over points in x-x' plane
for i=1:Ni
   %tic
   %     for j=1:Ni       
   %         %Find coordinates in kt-r plane
   %         kt=X(i,j);
   %         r=Xdash(i,j)/X(i,j);
   %         
   %         if r<=R(end,1) & r>=R(1,1); 
   %             %Find index of corrsponding to value of kt and r
   %             KTmap(i,j)=interp1(KT(1,:),1:size(KT,2),kt);
   %             Rmap(i,j)=interp1(R(:,1),1:size(R,1),r);
   % j        end
   %     end        
   %toc
   
   %tic
   kt=X(i,:);
   r=Xdash(i,:)./X(i,:);
   where=(r<=R(end,1) & r>=R(1,1));
   
   KTmap(i,where)=interp1(KT(1,:),1:size(KT,2),kt(where));
   Rmap(i,where)=interp1(R(:,1),1:size(R,1),r(where));
   
   %toc
end



%=============================================================
%=============================================================

function M=UETCmatrix(KTmap,Rmap,C,noDataValue,symmetry)

if nargin<5; symmetry=1; end
if nargin<4; noDataValue=1e-10; end
Ni=size(KTmap,1);
M=ones(Ni,Ni)*noDataValue;

%Loop over points in x-x' plane
for i=1:Ni
   for j=1:Ni
      
      if KTmap(i,j)~=0;
         %Find index of nearest lower value of KT and R
         ktIndex=floor(KTmap(i,j));
         rIndex=floor(Rmap(i,j));
         ktFrac=KTmap(i,j)-ktIndex;
         rFrac=Rmap(i,j)-rIndex;
         
         %Interpolate
         %M1 = linear Talyor expansion result from low kt, low r point
         %M2 = linear Talyor expansion result from high kt, high r point
         Maa=C(rIndex,ktIndex);
         Mab=C(rIndex,ktIndex+1);
         Mba=C(rIndex+1,ktIndex);
         Mbb=C(rIndex+1,ktIndex+1);
         M1=Maa+rFrac*(Mba-Maa)+ktFrac*(Mab-Maa);
         M2=Mbb-(1-rFrac)*(Mbb-Mab)-(1-ktFrac)*(Mbb-Mba);
         M(i,j)=(M1+M2)/2;

      end      
   end
end

%Put in symmetry condition (if required)
if symmetry==1
   for i=1:Ni
      for j=i:Ni
         if M(i,j)~=noDataValue & M(j,i)~=noDataValue
            Maverage=(M(i,j)+M(j,i))/2;
            M(i,j)=Maverage;
            M(j,i)=Maverage;
         elseif M(i,j)==noDataValue & M(j,i)~=noDataValue
            M(i,j)=M(j,i);
         elseif M(j,i)==noDataValue & M(i,j)~=noDataValue
            M(j,i)=M(i,j);
         end    
      end
   end
end

%=============================================================
%=============================================================

function M=UETCmatrixScalar(KTmap,Rmap,C11,C12,C21,C22)

noDataValue=1e-10;
Ni=size(KTmap,1);

%First get matrices ignoring coupling
M11=UETCmatrix(KTmap,Rmap,C11);
M12=UETCmatrix(KTmap,Rmap,C12,noDataValue,0);
M21=UETCmatrix(KTmap,Rmap,C21,noDataValue,0);
M22=UETCmatrix(KTmap,Rmap,C22);

%Insert symmetry for M12 and M21
for i=1:Ni
   for j=1:Ni
      if M12(i,j)~=noDataValue & M21(j,i)~=noDataValue
         Maverage=(M12(i,j)+M21(j,i))/2;
         M12(i,j)=Maverage;
         M21(j,i)=Maverage;
      elseif M12(i,j)==noDataValue & M21(j,i)~=noDataValue
         M12(i,j)=M21(j,i); 
      elseif M21(j,i)==noDataValue & M12(i,j)~=noDataValue
         M21(j,i)=M12(i,j);  
      end
   end
end

%Form full matrix
M=[M11 M12; M21 M22];

%=============================================================
%=============================================================

function output(file,kt,Evalue,Evector,scalar,Evaldim,vector)%JO Evaldim added

Ni=size(kt,2);

    j=1:Evaldim;
    k=1:Evaldim;
    [J,K]=meshgrid(j,k);
        
    
if Evaldim>512
   %disp(['Truncating eValues to first 512...'])
   %Evalue=Evalue(1:512);
%else 
   Evalue=Evalue(1:Evaldim);
end

if scalar==1
   [ktPhi eigenPhi]=eigenExtrapolate(kt,Evector(1:Ni,:),Ni,1);
   [ktPsi eigenPsi]=eigenExtrapolate(kt,Evector(Ni+1:end,:),Ni,1); 
   disp(['ktI  ' num2str(size(ktPhi))])
   disp(['vector  ' num2str(size(eigenPhi))])
   disp(['Evalue  ' num2str(size(Evalue))])
   M=[2*max(size(ktPhi)) ktPhi ktPsi; Evalue [eigenPhi; eigenPsi]'];
else
   [kt,Evector]=eigenExtrapolate(kt,Evector,Ni,0);
   M=[max(size(kt)) kt; Evalue Evector'];
end

tic
save(file,'-ascii','M');
t=toc;
disp(['Written: ' file ' in ' num2str(t) ' seconds']);

%=============================================================
%=============================================================

function [ktI,EvectorI]=eigenExtrapolate(kt,Evector,Evaldim,scalar) %Jo Evaldim added

tic;

N=size(Evector,2);

%if N>512
   disp(['Truncating to first' num2str(Evaldim) 'eigenvectors...'])
   Evector=Evector(:,1:Evaldim);
%end

%Extrapolate from kt(end) onwards
%ktExtra = [kt(end) + (kt(end)-kt(end-1))*[1:50]];
%EvectorExtra=repmat(Evector(end,:),50,1).*repmat([exp(-linspace(1,10,50))]',1,size(Evector,2));
%kt=[kt ktExtra];
%Evector=[Evector; EvectorExtra];

%Perform piecewise cubic Hermite interpolation

if scalar ==1
ktI=logspace(log10(kt(1)),log10(kt(end)),8*N);
else
ktI=logspace(log10(kt(1)),log10(kt(end)),4*N);
end

%figure
for i=1:min([size(Evector,2) Evaldim])
   EvectorI(:,i)=interp1([kt],[Evector(:,i)],ktI','pchip');
   % semilogx(kt,Evector(:,i),'k',ktI,EvectorI(:,i),'r'); hold on
end

%max(max(abs(Evector)))

t=toc;
disp(['Extrapolation and interpolation took: ' num2str(t) ' seconds'])

%=============================================================
%=============================================================

function D=UETCextrap(C,kt,r,ETC,ETC2)

%Hard-wired setting (should make an input???)
ktExtrap = 80; %Measured at reference time

disp(['Performing power-law extrapolation from kt = ' num2str(ktExtrap) ' (hardcoded)'])
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

%figure(15)
%clf
%loglog(kt,abs(C(1,:)),'k',kt,abs(D(1,:)),'b','LineWidth',2); hold on
%loglog(kt,abs(C(2,:)),'r',kt,abs(D(2,:)),'g'); hold on
%loglog(kt,abs(C(3,:)),'m',kt,abs(D(3,:)),'c'); hold on
%waitforbuttonpress

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