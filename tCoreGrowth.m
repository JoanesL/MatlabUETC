function [tCG rDiss] = tCoreGrowth(tEnd,tDiss,coreGrowthIndexA,coreGrowthIndexB,eraA,eraB,rCG)
%
% [tCG rCG] = tCoreGrowth(tEnd,tDiss,coreGrowthIndexA,coreGrowthIndexB,eraA,eraB,rCG)
%
% Calculates tCoreGrowth (tCG) for LAH, where tCoregrowth is the time 
% of switching between first core growth phase (A) and second (B).
% Also gives rCG  (width(tCG)/width)
%
%  tEnd  - end time of simulation
%  tDiss - end of dissipation (assumed same as tDiffusion) and start of
%          first core growth phase.
%  eraB  - power of t in scale factor a(t) = (t/tEnd)^eraB
%  eraA  - power of t in scale factor a(t) = a(tCG)*(t/tCG)^eraA
%  coreGrowthIndexB - width(t) = width*a(t)^coreGrowthIndexB
%  coreGrowthIndexA - width(t) = width(tCG)*[a(t)/a(tCG)]^coreGrowthIndexA
%  rCG   - width(tCG)/width
%

% if ~exist('rDiss','var')
%     rDiss = 1;
% end

tStart = tDiss-20; % Assumption!

pA = eraA*coreGrowthIndexA;
pB = eraB*coreGrowthIndexB;

%invP1 = 1/(pB-pA);
%invP2 = 1 - eraB*coreGrowthIndexB/(eraA*coreGrowthIndexA);

% tCG = (rDiss)^invP1*(tDiss/tEnd)^(1/invP2)*tEnd;
tCG = (rCG)^(1/pB)*tEnd;

%aCG = (tCG/tEnd)^eraB
%aDiss = aCG*(tDiss/tCG)^eraA

% rCG = aCG^coreGrowthIndexB;
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

plot(t(t>tStart),width(t>tStart),t(t<tCG),0.025*t(t<tCG))
legend('String width relative to end','\xi/10 estimate')

