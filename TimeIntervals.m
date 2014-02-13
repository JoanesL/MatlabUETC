function Intervals=TimeIntervals(q)

paperTimes=[0.049 0.150 0.308 0.536 0.852 1.207 1.609 2.069 2.6 3.219 3.95 4.828 5.901 7.243 8.967 11.27 14.49 19.31 27.36 43.46 91.74 156.1 478];
%EEEOOOO=paperTimes*123.9166
%paperx=1:1:22;

Intervals=logspace(log10(0.025),log10(paperTimes(end)),q);

Intervals=[0 Intervals];
end