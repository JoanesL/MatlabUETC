function PlotLike(path)

x=0.13:0.001:0.25;
y=0.1:0.001:0.2;

file=[path 'B2_logLike.txt'];

like=load(file);

datuk=like(:,3);

z=zeros(size(x,2),size(y,2));

l=1;

for i=1:size(x,2)
    for j=1:size(y,2)
        z(i,j)=datuk(l);
        l=l+1;
    end
end

contourf(x',y',z',500)