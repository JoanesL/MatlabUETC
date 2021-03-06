function PlotLike(path)

% x=0.13:0.001:0.25;
% y=0.1:0.001:0.2;
% 
% file=[path 'B2_logLike.txt'];

file_r=[path 'B2_3yr_camb_planck_withB_uK_20140314.txt'];
r_sp=load(file_r);

TT_sp=r_sp(:,2);

file_tex=[path 'SUM_001-128_00.dat'];
tex_sp=load(file_tex);

TT_tex=tex_sp(:,2);

T_cmb=2.726;
Factor=(T_cmb^2)/(2*pi);

TT_tex=Factor*TT_tex;

f10_factor = TT_tex(9) / TT_sp(9);

y=0:0.005:0.6;
y=y*f10_factor;
x=0:0.01:0.33;

file=[path 'B2_logLike_tex.txt'];

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

% temp=max(z);
% Zmax=max(temp)
% 
% z=z/Zmax;

z=z*(-2); %\chi^2

temp=min(z);
zmin=min(temp)

z=z-zmin;


% temp=max(z);
% Zmax=max(temp);
% 
% z=z/Zmax;
% 
contour(x',y',z',[2.3 6.17])
% 
% sumZ2= sum(z,2)
% sumZ1= sum(z,1)
% 
% figure()
% plot(x,sumZ2)
% figure()
% plot(y,sumZ1)

%contourf(x',y',z',10)
