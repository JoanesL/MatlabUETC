function ploterrorbars(x,y,e_m,e_p,string)

if nargin==4
   string=e_p;
   e_p=e_m;
end

next_plot=get(gca,'NextPlot');
hold on;

for i=1:max(size(x))
   plot(x(i)*ones(1,3),y(i)*ones(1,3)+[-e_m(i) 0 e_p(i)],['-' string],'LineWidth',1)
end
set(gca,'NextPlot',next_plot);
