%% Caculation of diversity index
function idx = d_index(P,LB,UB)
    NP = length(P(:,1));
    Bound = UB - LB; Id=Bound~=0;
    xc = mean(P);
    idx = 0;
    for i=1:NP
        x = P(i,:);
        idx = idx + (sum(((xc(Id)-x(Id))./Bound(Id)).^2))^0.5;
    end
    idx = idx/NP;
end