%% Sorting of constrained solution
function [S,ID]=consort(F,C,k,e)
  S = zeros(1,length(F));
  for i=1:length(F)
      if C(i)<=e, S(i)=F(i)/(F(i)+k);
      else S(i)=1+C(i);
      end
  end
  [S,ID]=sort(S);
end