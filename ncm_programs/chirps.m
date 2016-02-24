% Bird chirps

load chirp
j = 1:1560;
for k = 1:8
   sound(y((k-1)*1560+j),Fs);
   subplot(4,2,k)
   p = abs(fft(y((k-1)*1560+j)));
   plot(p(1:780))
   title(num2str(k))
   axis([390 780 0 40])
   set(gca,'xtick',[390 780],'ytick',[])
   drawnow
   pause
end
