function [lambda charge] = ccEnd(tStart,tDiss,tCG,tEnd_s,tEnd,coreGrowthIndexA,coreGrowthIndexB,eraA,eraB)
%
% [lambda charge] = ccEnd(tStart,tDiss,tCG,tEnd_s,coreGrowthIndexA,coreGrowthIndexB,eraA,eraB)
%
% Calculates values of scalar and gauge coupling constants at tEnd_s, the
% end time of a short simulation intended to test initial conditions.
%
%  tStart  - start time of simulation
%  tDiss - end of dissipation (assumed same as tDiffusion) and start of
%          first core growth phase.
%  tCG - tCoreGrowth
%  tEnd_s   - end time of small simulation
%  tEnd  - end time of simulation
%  eraB  - power of t in scale factor a(t) = (t/tEnd)^eraB
%  eraA  - power of t in scale factor a(t) = a(tCG)*(t/tCG)^eraA
%  coreGrowthIndexB - width(t) = width*a(t)^coreGrowthIndexB
%  coreGrowthIndexA - width(t) = width(tCG)*[a(t)/a(tCG)]^coreGrowthIndexA
%  rCG   - width(tCG)/width
%


pA = eraA*coreGrowthIndexA;
pB = eraB*coreGrowthIndexB;

aCG = (tCG/tEnd)^eraB;
aDiss = aCG*(tDiss/tCG)^eraA;

rCG = aCG^coreGrowthIndexB;
rDiss = rCG*(tDiss/tCG)^pA;

numt = 500;

t = linspace(0,tEnd,numt);

width = zeros(size(t));

for n = 1:numt
    if t(n) <= tDiss
        width(n) = rDiss;
    end
    if t(n) > tDiss && t(n) <= tCG
        width(n) = rCG*(t(n)/tCG)^pA;
    end
    if t(n) > tCG
        width(n) = (t(n)/tEnd)^pB;
    end
end

[t width] = rCoreGrowth(tStart,tDiss,tCG,tEnd,coreGrowthIndexA,coreGrowthIndexB,eraA,eraB);

wEnd_s = interp1(t,width,tEnd_s);

lambda = 2/wEnd_s^2;
charge = 1/wEnd_s;
plot(t(tStart <t & t<tEnd_s),width(tStart <t & t<tEnd_s),t(t<tCG),0.025*t(t<tCG))
legend('String width','\xi/10 estimate')

