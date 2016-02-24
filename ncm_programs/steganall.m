% STEGANALL  Steganography.  Here are all images hidden in the
% default CDATA for the IMAGE command.

p = [1  6 11 16 17 18 19 20 25 30 35 36 40 44 48 52];
q = [5 10 15 16 17 18 19 24 29 34 35 39 43 47 51 52];

clf
image
imageh = get(gca,'child');
cdata = get(imageh,'cdata')/32;
clf
shg
colormap(gray(32))
for k = 1:16
   subplot(4,4,k)
   m = p(k);
   n = q(k);
   imagesc(mod(floor(2^n*cdata),2^(n-m+1))+1)
   title([int2str(m) ':' int2str(n)])
   axis image
   axis ij
   axis off
end
