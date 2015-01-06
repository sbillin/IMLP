function [v] = EddingtonEps(i,j,k)
  v = (j-i)*(k-j)*(k-i)/2;
end
