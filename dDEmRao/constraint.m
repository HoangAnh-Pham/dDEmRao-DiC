%% Consntraint violation 
function cv = constraint(fcons,x,delta)
if isempty(fcons)==1, cv = 0;
else
[C, Ceq] = feval(fcons,x);
Co = [C,Ceq];
q = length(C);
m = q+length(Ceq);
G = zeros(1,m);
for i=1:m
    if i<=q, G(i) = max(0,Co(i));
    else
        G(i) = max(0,abs(Co(i))-delta);
        %G(i) = max([0,Co(i),-Co(i)]);
    end
end
%cv = sum(G);
%cv = max(G);
cv = max(Co);
end
end