% Read the text of Lincoln's Gettysburg address with

   fp = fopen('gettysburg.txt');
   G = char(fread(fp))'
   fclose(fp);

% How many characters are in the text?
% Use the 'unique' function to find the unique characters in this text? 
% What punctuation characters, and how many of each, are there?

   nchar = length(G)
   uniq = unique(G)
   nuniq = length(uniq)
   nblank = sum(G==' ')
   nperiod = sum(G=='.')
   ncomma = sum(G==',')
   ndash = sum(G=='-')

% Remove the punction and convert the text to one case.

   G = upper(G);
   G(G<'A') = [];
   G(G>'Z') = [];

% Use the 'histc' function to find the most frequent letter
% and the missing letters.

   G = double(G);
   b = double('A':'Z');
   N = histc(G,b);
   most = char(b(find(N==max(N))))
   missing = char(b(find(N==0)))

% Use the 'bar' function as described in 'help histc'
% to plot a histogram of the letter frequencies.
% Use get(gca,'xtick') and get(gca,'xticklabel') to see
% how  the x-axis of the histogram is labeled.  Then use
% set(gca,'xtick',...,'xticklabel',...)
% to label the x-axis with the letters in the text.

   bar(b,N,'histc')
   set(gca,'xlim',[min(b)-2 max(b)+2])
   xtick = get(gca,'xtick');
   xticklabel = get(gca,'xticklabel');
   x = b;
   x(find(N==0)) = [];
   set(gca,'xtick',x+.5,'xticklabel',char(x'))
   title('Gettysburg Address')
