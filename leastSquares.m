function p=leastSqu(x,y)

N = max(size(x));
Sx = sum(x);
Sy = sum(y);
Sxx = sum(x.*x);
Sxy = sum(x.*y);

D = N*Sxx - Sx^2;
m = (N*Sxy - Sx*Sy) / D;
c = (Sxx*Sy - Sx*Sxy) / D;

p = [m c];   %polyfit format answer