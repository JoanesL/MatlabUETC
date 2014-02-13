
%FUNCTION PLOTC(X,Y,V,'MARKER') plots the values of v colour coded
% at the positions specified by x and y, and v (z-axis) in a 3-D axis
% system. A colourbar is added on the right side of the figure.
%
% The colorbar strectches from the minimum value of v to its
% maximum in 9 steps (10 values).
%
% The last argument is optional to define the marker being used. The
% default is a point. To use a different marker (such as circles, ...) send
% its symbol to the function (which must be enclosed in '; see example).
%
% The plot is actually a 3D plot but the orientation of the axis is set
% such that it appears to be a plane 2D plot. However, you can toggle
% between 2D and 3D view either by using the command 'view(3)' (for 3D
% view) or 'view(2)' (for 2D), or by interactively rotating the axis
% system.
%
% Example:
% Define three vectors
%    x=1:10;y=1:10;p=randn(10,1);
%    plotc(x,y,p)
%
%    x=randn(100,1);
%    y=2*x+randn(100,1);
%    p=randn(100,1);
%    plotc(x,y,p,'d')
%    view(3)
%
% Uli Theune, University of Alberta, 2004
% modified by Stephanie Contardo, British OCeanographic Data Centre, 2006
% input scalar=1 if dotproduct of scalar pert, otherwise scalar=0

function plotc(x,y,v,Evalmin,Nvalue,markersize,interp,marker)
if nargin==0; 
  help plotc
  return
end

delete(gca)
if nargin <8
    marker='o';
    markersize=20;
end
%colormap(lines(128))
map=colormap;
%Jo z-ren balio hauetarako marraztuko ditu
if interp==1
    miv=-50;
    mav=50;
else
miv=-2.05;%min(miv1);
mav=2.05;%max(mav1);
end
clrstep = (mav-miv)/size(map,1); 
col=linspace(miv,mav,10);
% Plot the pointshot)%
hold on

%for i=Evalmin:Nvalue
for i=1:(Nvalue+1-Evalmin)
for nc=1:size(map,1)
     
        
    iv = find(v(i,1:end)>miv+(nc-1)*clrstep & v(i,1:end)<=miv+nc*clrstep) ;
  
    %plot3(x(1:end,i),y(1:end,1),v(i,1:end),marker,'color',map(v(i,i),1:end),'markerfacecolor',map(v(i,i),1:end))
    plot3(x(iv,i),y(iv,1),v(i,iv),marker,'color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',markersize)
    set(gca,'Ylim',[Evalmin Nvalue],'YTick',Evalmin:Nvalue)
    set(gca,'Xlim',[Evalmin Nvalue],'XTick',Evalmin:Nvalue)
end
end
      
hold on 

% Re-format the colorbar
h=colorbar;

set(h,'ylim',[1 length(map)]);
yal=linspace(1,length(map),15);
set(h,'ytick',yal);
 %Create the yticklabels
ytl=linspace(miv,mav,15);
s=char(clrstep,15);
for i=1:15
    if min(abs(ytl)) >= 0.001
        B=sprintf('%-4.3f',ytl(i));
    else
        B=sprintf('%-3.1E',ytl(i));
    end
    s(i,1:length(B))=B;
end
set(h,'yticklabel',s);
grid on
view(2)
