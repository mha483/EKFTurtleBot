clear, clc

r = [ 1 4 6 7;
      2 4 6 1;
      3 0 4 7;
      4 8 6 9;
      5 5 0 8;
      6 7 8 9;
      7 1 1 1];

rowsToDel = max(r==0,[],2);
r(rowsToDel,:) = [];

