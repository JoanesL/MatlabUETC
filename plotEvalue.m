function plotEvalue()

S=load(['../UETCdata/UETCdata2013/Merged/Baseline/UETCInterpolation11/Interpolated/UETCeigenVector_01_11.dat']);

for i=1:100
v.Evalue(i) = S(i+1,1);
end

figure()
plot(1:100, v.Evalue, 'b', 1:100, abs(v.Evalue),'k')