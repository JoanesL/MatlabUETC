function [kt,r,C1,C2]=UETCplot(Cname,id,run,tRef,tOffset,xiscaling,inPath,parent)
%UETC plotting function for UETC.hpp data 2006-2013
%
% Usage: [kt,r,C1,C2]=UETCplot(Cname,id,run,tRef,tOffset,inPath,parent)
%
%  Cname = UETC name, eg. scalar11 or vector
%          (add * at start to plot coherence plot, ie. /ETC)
%     id = ID string between 'UETCscalar11_' and before '.dat'
%          eg. to load statsFile_6L01 ID is '6L%2' with run=1
%    run = realizations(s) to include
%   tRef = UETC reference time
%tOffset = time when xi=0, if '*' get from statsFile (Lag & tRef<t<2*tRef)
%xiscaling = 1-> xi scaling, toffset=0
%            0-> toffset rescaling
%Optional parameters are:
%
%   inpath = path to file, including final '/'
%          (if ommited or '' gets path from gpath global variable)
% parent = destination axes for plot
%
% 13.10.04 MBH tOffset now performed for each run rather than on average
%              Can now handle data files without runID
%
% TO DO: allow xi scaling (as with ETCplot/load) in UETCtimeOffset


if nargin==0; 
  help UETCplot
  return
end

global gpath

if ~exist('inPath','var'); inPath=''; end 

if prod(size(inPath))>0; 
  path=inPath; 
else
  if prod(size(gpath))>0
    path=gpath;
  else
    disp(['Please set gpath global variable to default path'...
	  ' or specify path in fucntion call'])
    return
  end
end

%Get coherence or normal UETC plot
if Cname(1)=='*'
  Cname=Cname(2:end);
  plotCoherence = 1;
  disp('Coherence plot, ie. UETC/ETC requested...'); 
else
  plotCoherence = 0;
end

%if strcmp(tOffset,'xiscaling')==1
%    xiscaling=1;
%end

%Take average of tOffsets if multiple = UETC plots for a common mean
%tOffSet = mean(tOffSet);

%Prepare axes
%if exist('parent','var')~=1; clf; else axes(parent); end  

%LOAD DATA
[kt,r,C1]=UETCload(path,Cname,id,run,tRef,tOffset,xiscaling);
% if exist('tOffSet','var')==1
%     disp(['Initial max ratio: ' num2str(max(r))])
%     [kt,r,C1]=UETCtimeOffSet(Cname,kt,r,C1,tRef,tOffSet); %Apply time offset correction
%     disp(['Post-offset max ratio: ' num2str(max(r))])
% end

if strfind(Cname,'12')>0
    [kt,r,C2]=UETCload(path,'scalar21',id,run,tRef,tOffset,xiscaling);
%     if exist('tOffSet','var')==1
%     [kt,r,C2]=UETCtimeOffSet('scalar21',kt,r,C2,tRef,tOffSet);  %Apply time offset correction
%     end
else
    C2=C1;
end

%Divide by ETC if coherence plot
if plotCoherence==1
   for i=2:size(r,1)
     C1(i,:) = C1(i,:) ./ C1(1,:);
     C2(i,:) = C2(i,:) ./ C2(1,:);
   end
   C1(1,:)=1;
   C2(1,:)=1; 
end

if strfind(Cname,'ten')>0
disp('Note that tensor data is not reduced by factor of 2')
disp('in this version of UETCplot');
end

%Limit kt (if desired)

%which=find( kt>49 & kt<=55);

%C1=C1(:,kt<=20);
%C1=C1(:,which);
%C2=C2(:,kt<=20);
%C2=C2(:,which);
%kt=kt(kt<=20);
%kt=kt(which);

    %zein1 = find( kXi1>=fitlimitkxi(1) & kXi1<=fitlimitkxi(2) );


%Limit r (if desired)
%C1=C1(r<1.3333333,:);
%C2=C2(r<1.3333333,:);
%r=r(r<1.33333333);

%CALCULATE R AND Z MATRICES

for i=1:size(kt,2)
    for j=1:size(r,1)
        R1(j,i)=r(j);
        R2(j,i)=1/r(j);
        Z1(j,i)=kt(i)*sqrt(r(j));
        Z2(j,i)=kt(i)*sqrt(r(j));
        X(j,i)=kt(i);
        Y(j,i)=kt(i)*r(j);
        Rprueba(j,i)=kt(i);
    end
end

%EXTRAPOLATE TO SUPERHORIZON SCALES
%if size(strfind(path,'gO4'),2)==0
%   for j=1:size(r,1)
%      zExtra=logspace(log10(1),log10(Z1(j,1)),2);
%      zExtra=zExtra(1:end-1);
%      Z1extra(j,:)=[zExtra Z1(j,:)];
%      R1extra(j,:)=[repmat(R1(j,1),size(zExtra)) R1(j,:)];
%      C1extra(j,:)=[repmat(mean(C1(j,2:3)),size(zExtra)) mean(C1(j,2:3)) C1(j,2:end)];
%      Z2extra(j,:)=[zExtra Z2(j,:)];
%      R2extra(j,:)=[repmat(R2(j,1),size(zExtra)) R2(j,:)];
%      C2extra(j,:)=[repmat(mean(C2(j,2:3)),size(zExtra)) mean(C2(j,2:3)) C2(j,2:end)];
%   end
%   Z1=Z1extra; R1=R1extra; C1=C1extra;
%   Z2=Z2extra; R2=R2extra; C2=C2extra;
%else
%   for j=1:size(r,1)
%      zExtra=logspace(log10(min(min(Z1))),log10(Z1(j,1)),2);
%      zExtra=zExtra(1:end-1);
%      Z1extra(j,:)=[zExtra Z1(j,:)];
%      R1extra(j,:)=[repmat(R1(j,1),size(zExtra)) R1(j,:)];
%      C1extra(j,:)=[repmat(mean(C1(j,2:3)),size(zExtra)) mean(C1(j,2:3)) C1(j,2:end)];
%      Z2extra(j,:)=[zExtra Z2(j,:)];
%      R2extra(j,:)=[repmat(R2(j,1),size(zExtra)) R2(j,:)];
%     C2extra(j,:)=[repmat(mean(C2(j,2:3)),size(zExtra)) mean(C2(j,2:3)) C2(j,2:end)];
%   end
%   Z1=Z1extra; R1=R1extra; C1=C1extra;
%   Z2=Z2extra; R2=R2extra; C2=C2extra;   
%end

figure;
colour=colourMap(C1);
%colour(Z1>100)=NaN;

surf(R1,Z1,abs(C1),colour,'FaceLighting','phong','FaceColor','interp','EdgeColor','flat');

%surf(R1,Z1,abs(C1),colour,'FaceLighting','phong','FaceColor','interp','EdgeColor','flat');
%surf(R1,Z1,abs(C1),'FaceLighting','phong','FaceColor','interp','EdgeColor','flat')
%colormap(flipud(gray));
colour=colourMap(C2);
%colour(Z2>100)=NaN;

hold on
surf(R2,Z2,abs(C2),colour,'FaceLighting','phong','FaceColor','interp','EdgeColor','flat');

%surf(R2,Z2,abs(C2),colour,'FaceLighting','phong','FaceColor','interp','EdgeColor','flat');
%surf(R2,Z2,abs(C2),'FaceLighting','phong','FaceColor','interp','EdgeColor','flat')
%colormap(flipud(gray));
%shading interp
view([-141 36])
%view([90 0])
%camlight right
%camlight left
camlight headlight
set(gca,'YScale','log','XScale','log')
set(gca,'ZScale','linear')
xlabel('t`/t')
ylabel('K (tt`)^{1/2}')
if(xiscaling == 1)
    if strcmp(Cname,'scalar11')==1
        zlabel('((\xi/t)(\xi`/t`))^{1/2} |C^s_{11}|')
    elseif strcmp(Cname,'scalar12')==1
        zlabel('((\xi/t)(\xi`/t`))^{1/2} |C^s_{12}|')    
    elseif strcmp(Cname,'scalar21')==1
        zlabel('((\xi/t)(\xi`/t`))^{1/2} |C^s_{21}|')
    elseif strcmp(Cname,'scalar22')==1
        zlabel('((\xi/t)(\xi`/t`))^{1/2} |C^s_{22}|')
    elseif strcmp(Cname,'vector')==1
        zlabel('((\xi/t)(\xi`/t`))^{1/2} |C^s_{vv}|')
    elseif strcmp(Cname,'tensor')==1
        zlabel('((\xi/t)(\xi`/t`))^{1/2} |C^s_{tt}|')
    end
elseif Cname(1)=='*'
    if strcmp(Cname,'scalar11')==1
        zlabel('|D^s_{11}|')
    elseif strcmp(Cname,'scalar12')==1
        zlabel('|D^s_{12}|')   
    elseif strcmp(Cname,'scalar21')==1
        zlabel('|D^s_{21}|')
    elseif strcmp(Cname,'scalar22')==1
        zlabel('|D^s_{22}|')
    elseif strcmp(Cname,'vector')==1
        zlabel('|D^s_{vv}|')
    elseif strcmp(Cname,'tensor')==1
        zlabel('|D^s_{tt}|')
    end
else
    if strcmp(Cname,'scalar11')==1
        zlabel('|C^s_{11}|')
    elseif strcmp(Cname,'scalar12')==1
        zlabel('|C^s_{12}|')   
    elseif strcmp(Cname,'scalar21')==1
        zlabel('|C^s_{21}|')
    elseif strcmp(Cname,'scalar22')==1
        zlabel('|C^s_{22}|')
    elseif strcmp(Cname,'vector')==1
        zlabel('|C^s_{vv}|')
    elseif strcmp(Cname,'tensor')==1
        zlabel('|C^s_{tt}|')
    end
end
    
set(gcf,'Color',[1 1 1]);
set(gca,'Color',[0.8 0.8 0.8]);

set(gca,'LineWidth',2);
set(gca,'FontSize',14)

rMax=6;
ktMin=kt(1);
ktMax=max(Z1);
ktMax=max(ktMax);
CMax=max(C1);
CMax=max(CMax)

set(gca,'XLim',[1/rMax rMax],'XTick',[0.125 0.25 0.5 1 2 4 8])
set(gca,'YLim',[ktMin ktMax])
%set(gca,'ZLim',[0 CMax])

if strfind(path,'gO4')>0
	set(gca,'XLim',[1/7 7],'XTick',[0.25 0.5 1 2 4])
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

function D=UETCextrap(C,kt,r,ETC,ETC2)

%Hard-wired setting (should make an input???)
ktExtrap = 200; %Measured at reference time

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
  for j=2:size(ETC,1)
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

%C=1-C;