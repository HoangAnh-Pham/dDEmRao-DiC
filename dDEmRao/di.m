function d = di(X,Y,LB,UB)
  Bound = UB*(1+0.000001)-LB;
  d = sqrt(sum(((X-Y)./Bound).^ 2));
end
