%% Simple bound constraints
function sol = boundConst(sol,x,Lb,Ub)
  ns_tmp=sol;
  % Apply the lower bound
  I=ns_tmp<Lb; 
  ns_tmp(I)=(x(I) + Lb(I)) / 2; 
  %ns_tmp(I)=Lb(I)+(Lb(I)-ns_tmp(I))-fix((Lb(I)-ns_tmp(I))/(Ub(I)-Lb(I)))*(Ub(I)-Lb(I)); 
  % Apply the upper bounds 
  J=ns_tmp>Ub; 
  ns_tmp(J)=(x(J) + Ub(J)) / 2;
  %ns_tmp(J)=Ub(J)-(ns_tmp(J)-Ub(J))+fix((ns_tmp(J)-Ub(J))/(Ub(J)-Lb(J)))*(Ub(J)-Lb(J)); 
  % Update this new move 
  sol=ns_tmp;
end
