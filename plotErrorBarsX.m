function ploterrorbars(x,y,e_m,e_p,string)

if nargin==4
   string=e_p;
   e_p=e_m;
end

next_plot=get(gca,'NextPlot');
hold on;

for i=1:max(size(x))
   plot(x(i)*ones(1,2)+[-e_m(i) e_p(i)],y(i)*ones(1,2),['-' string])
end
set(gca,'NextPlot',next_plot);