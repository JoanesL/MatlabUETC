%Usage: UETCeigenInterp(idCell,runCell,tRef,tOffSet,pathCell,outPath,eraCell,number,bootstrap,Ni,rMax,plateauMode,highKtMode)
%
%   idCell = ID string between 'UETCscalar11_' and before '.dat'
%               e.g., Merged Baseline: {_mat_%2,_rad_%2}
%    runCell = realizations(s) to include. e.g {1,1}
%     tRef = UETC reference time. 
%  tOffset = time when xi=0, if '*' get from statsFile (Lag & tRef<t<2*tRef)
%            (independent offsets for each realization is used, 
%             but a common one (ie. single number) can be specified too
%
% pathCell = path to input UETCs files, including final '/'
%            (if ommited or '' gets path from gpath global variable)
%  outPath = path (relative to inPath) for outputting UETCeigen files
%  eraCell = either 'Mat' or 'Rad' (simply for output file names)
%   number = ID number in output filenames, default 1
%bootstrap = 0 means no, use average 
%            1 means use average of random 
%            realizations combinations (without removal)
%       Ni = (interpolated) matrix size (default 512)
%     rMax = zero data above rMax
%plateauMode = 0 means use lowest kt point (point 1)
%              1 means use mean of points 2 and 3 (and replace 1)
%              (default is 1)
%highKtMode = 0 means use raw data for high Kt
%             1 means extrapolate/correct as power-law at high Kt
%example of Merged Baseline:
%UETCeigenInterp({'_mat_%2','_rad_%2'},{1,1},[600 450],0,{pathMerged_Baseline,pathMerged_Baseline},'UETCInterpolation22/',{'Mat','Rad'},1,0,2048,10,4,0)

function UETCeigenInterp(idCell,runCell,tRef,tOffSet,pathCell,outPath,eraCell,number,bootstrap,Ni,q,rMax,plateauMode,highKtMode,negativeOrder)

if nargin==0; 
  help UETCeigenInterp
  return
end

global gpath

if ~exist('inPath','var'); inPath=''; end 

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
% Convert to cell array even for one element
for n = 1:nPaths
    tOff{n} = tOffSet;
end  

if ~exist('scale','var')
    scale = 'log';
end


if ~exist('outPath','var'); outPath=''; end
if ~exist('bootstrap','var'); bootstrap=0; end
if ~exist('Ni','var'); Ni=512; end
if ~exist('rMax','var'); rMax=99; end
if rMax==0; rMax=99; end
if ~exist('plateauMode','var'); plateauMode=1; end
if ~exist('highKtMode','var'); highKtMode=1; end

%Make output ID string
if prod(size(outPath))>0
    for i=1:numel(eraCell)
        era{n}=[eraCell{n} '_' num2str(number,'%2.2d')];
    end
end
  
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

%Get interpolation coordinates (x=kt)
disp(['Matrix size: ' num2str(Ni) 'x' num2str(Ni)])

%Joanes Lizarraga (28/11/13)
%Hardcoded for Merged UETCs
x=linspace(0.1,2350,Ni);

disp(['Matrix kt from:' num2str(min(x)) ' to ' num2str(max(x)) ' (hardcoded)'])

[X,Xdash]=meshgrid(x,x);

for n=1:numel(pathCell)
    runs=runCell{n};
    tOffSet=tOff{n};
    %Expand tOffSet to vector if only one given for multiple runs
    if prod(size(tOffSet))==1; tOffSet=tOffSet*ones(size(runs)); end

    for i=1:prod(size(runs))
        %LOAD CORRELATORS
        %If Merged UETCs DON'T USE xiScaling, they are already correctly
        %scaled
        [kt,r,s.C11]=UETCload(pathCell{n},'scalar11',idCell{n},runs(i),tRef(n),tOffSet(i),0); 
        [kt,r,s.C12]=UETCload(pathCell{n},'scalar12',idCell{n},runs(i),tRef(n),tOffSet(i),0);
        [kt,r,s.C21]=UETCload(pathCell{n},'scalar21',idCell{n},runs(i),tRef(n),tOffSet(i),0);
        [kt,r,s.C22]=UETCload(pathCell{n},'scalar22',idCell{n},runs(i),tRef(n),tOffSet(i),0);
        [kt,r,v.C]=UETCload(pathCell{n},'vector',idCell{n},runs(i),tRef(n),tOffSet(i),0);
        [kt,r,t.C]=UETCload(pathCell{n},'tensor',idCell{n},runs(i),tRef(n),tOffSet(i),0);
   
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
            [temp,temp,ETC1]=ETCload(pathCell{n},'scalar11',idCell{n},runs(i),tRef(n),tOffSet(i),0,999);
            s.C11=UETCextrap(s.C11,kt,r,ETC1);

            [temp,temp,ETC2]=ETCload(pathCell{n},'scalar22',idCell{n},runs(i),tRef(n),tOffSet(i),0,999);
            s.C22=UETCextrap(s.C22,kt,r,ETC2);

            s.C12=UETCextrap(s.C12,kt,r,ETC1,ETC2);
            s.C21=UETCextrap(s.C21,kt,r,ETC2,ETC1);

            [temp,temp,ETC1]=ETCload(pathCell{n},'vector',idCell{n},runs(i),tRef(n),tOffSet(i),0,999);
            v.C=UETCextrap(v.C,kt,r,ETC1);

            [temp,temp,ETC1]=ETCload(pathCell{n},'tensor',idCell{n},runs(i),tRef(n),tOffSet(i),0,999);
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
            s_M=UETCmatrixScalar(KTmap,Rmap,s.C11,s.C12,s.C21,s.C22);
            v_M=UETCmatrix(KTmap,Rmap,v.C);
            t_M=UETCmatrix(KTmap,Rmap,t.C);
        else
            s_M=s_M+UETCmatrixScalar(KTmap,Rmap,s.C11,s.C12,s.C21,s.C22);
            v_M=v_M+UETCmatrix(KTmap,Rmap,v.C);
            t_M=t_M+UETCmatrix(KTmap,Rmap,t.C) ;       
        end    
    end

    %AVERAGE MATRICES
    disp('AVERAGING MATRICES...')
    s_M=s_M/prod(size(runs));
    v_M=v_M/prod(size(runs));
    t_M=t_M/prod(size(runs));
    
    Cscalar{n}=s_M;
    Cvector{n}=v_M;
    Ctensor{n}=t_M;
end

%Hardcoded Dani-Martin paper
t_eq=123.916;
t_frac=[0 0.049 0.254 0.852 1.609 2.600 3.950 5.901 8.967 14.49 27.36 91.74 478 1000.0];
%t_frac_tarte=linspace(8.9665,8.967,50);
%t_frac_tarte2=linspace(8.966, t_frac_tarte(2),25);
%t_frac3=linspace(t_frac_tarte2(1),t_frac_tarte2(3),25);
%t_frac=linspace(8.966,t_frac3(2),50);
%t_frac=logspace(0,478,33);
%t_frac=[0 0.049 0.150 0.308 0.536 0.852 1.207 1.609 2.069 2.6 3.219 3.95 4.828 5.901 7.243 8.967 11.27 14.49 19.31 27.36 43.46 91.74 156.1 478 99999];
%t_frac=TimeIntervals(q);

%Save Time Intervals for CMBeasy
t_int=t_frac*t_eq;
filetimeInt=[pathCell{1} outPath 'UETCeigentimeIntervals_' num2str(number,'%2.2d') '.dat']
save(filetimeInt,'-ascii','t_int')

figure
a=multiPlot([2 1]);
for i=1:size(t_frac,2)
    
    a(i)=(((2^(1/2) - 1) * t_frac(i) + 1)^2 - 1);
    %f_neil(l)=(1+a(l))^(-1);
    
    disp(['Time interval ' num2str(i) ])
    
    if i==1
        f_t(i)=1;
    elseif i==size(t_frac,2)+1
        f_t(i)=0;
    else
        %f_t(i)=power((1+0.25*t_frac(i)),-1)
        %f_neil=(e_neil)^2
        f_t(i)=(1+a(i))^(-2);
    end
% for i=1:size(t_frac,2)-1
%     
%     disp(['Time interval ' num2str(i) ])
%     
%     if i==size(t_frac,2)
%         f_t(i)=0
%     elseif i==1
%         t_teq(i)=0;
%     else
%         t_teq(i)=(t_frac(i));%+t_frac(i+1))/2;
%     end
%     
%         f_t(i)=power((1+0.25*t_teq(i)),-1)
% %     end
    
    %n=2 Radiation, n=1 Matter
    s_M=f_t(i)*Cscalar{2} + (1-f_t(i))*Cscalar{1};
    v_M=f_t(i)*Cvector{2} + (1-f_t(i))*Cvector{1};
    t_M=f_t(i)*Ctensor{2} + (1-f_t(i))*Ctensor{1};
    
    %CALCULATE EIGENVALUES AND VECTORS
    disp('CALCULATING EIGENVALUES AND EIGENVECTORS...')
    [s_Evector,s_Evalue]=eigen(s_M,negativeOrder);
    [v_Evector,v_Evalue]=eigen(v_M,negativeOrder);
    [t_Evector,t_Evalue]=eigen(t_M,negativeOrder);
    

    %Evolution of first 3 evectors
%     axes(a(1))
%     semilogx(x, s_Evector(1:Ni,1),'r', x, s_Evector(1:Ni,2),'g', x, s_Evector(1:Ni,3),'b'); hold on 
%     semilogx(x, s_Evector(Ni+1:end,1),'--r', x, s_Evector(Ni+1:end,2),'--g', x, s_Evector(Ni+1:end,3),'--b'); hold on
%     xlabel('index')
%     ylabel('Eigenvector value')
%     axis tight
%     
%     axes(a(2))
%     plot(1:Ni, v_Evalue,'b', 1:Ni,abs(v_Evalue),'r', [1 Ni],[0 0],'k'); hold on;
%     xlabel('Eigenvalue number')
%     ylabel('Eigenvalue')
%     axis tight
%     set(gca,'Xlim',[1 50])
    
%     figure()
%     plot(1:Ni, v_Evalue,'b', 1:Ni,abs(v_Evalue),'r', [1 Ni],[0 0],'k'); hold on;
%     xlabel('Eigenvalue number')
%     ylabel('Eigenvalue')
%     set(gca,'Xlim',[1 50])
%     titleplot=['e-Values: ' num2str(i) ];
%     title(titleplot ,'FontWeight','bold')    
    

    %Save eigenvecotors: usual format
    fileeigensca=[pathCell{1} outPath 'UETCeigenScalarBefore_' num2str(number,'%2.2d') '_' num2str(i,'%2.2d') '.dat']
    fileeigenvec=[pathCell{1} outPath 'UETCeigenVectorBefore_' num2str(number,'%2.2d') '_' num2str(i,'%2.2d') '.dat'];
    fileeigenten=[pathCell{1} outPath 'UETCeigenTensorBefore_' num2str(number,'%2.2d') '_' num2str(i,'%2.2d') '.dat'];

    A=[2*max(size(x)) x x; s_Evalue (s_Evector)'];
    B=[max(size(x)) x; v_Evalue (v_Evector)'];
    C=[max(size(x)) x; t_Evalue (t_Evector)'];

    save(fileeigensca,'-ascii','A');
    save(fileeigenvec,'-ascii','B');
    save(fileeigenten,'-ascii','C');
    disp('Saving eigenvalues and eigenvectors...')

end



%=============================================================
%=============================================================

function [vector,value]=eigen(C,negativeOrder)

N=size(C,1);

%Use eig to find eigenvalues and vectors
tic
[vector,value]=eig(C);
t=toc;
disp(['Eigenvectors determined in ' num2str(t) 's'])
value=diag(value);

%Sort in descending eigenvalue (magnitude... if desired)
if (negativeOrder==1)
    [dummy,index]=sort(-abs(value));
else
    [dummy,index]=sort(-value);
end
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
ktI=logspace(log10(0.4),log10(2450),8*N);
else
ktI=logspace(log10(0.4),log10(2450),4*N);
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
  for j=1:max(size(r))
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

