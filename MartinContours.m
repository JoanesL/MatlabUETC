function MartinContours()
nc=2; % number of contours
ct = zeros(nc);
ct(2)=0.683;
ct(1)=0.954;
% ... get contours ...
norm = sum(sum(bin));
for j = 1:nc
   try_t = max(max(bin));
   try_b = 0;
   try_sum = 0;
   try_sav = -1;
   while (abs(try_sav-try_sum)~=0)
       try_sav = try_sum;
       mask = find(bin<(try_b+try_t)/2);
      try_sum = sum(bin(mask));
      if (try_sum> (1-ct(j))*norm)
         try_t = (try_b+try_t)/2;
      else
         try_b = (try_b+try_t)/2;
      end
  end % while loop
   cl(j) = (try_b+try_t)/2;
end
% .... generate axes ...
ax1=zeros(nb);
ax2=zeros(nb);
for ix1=1:nb
   ax1(ix1,:) = minp(1)+(ix1-0.5)*width(1);
   ax2(:,ix1) = minp(2)+(ix1-0.5)*width(2);
end
%figure(2)
contour(ax1,ax2,bin,cl,'k','LineWidth',2)