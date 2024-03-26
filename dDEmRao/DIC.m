%% Distance Comparison Method
function eval=DIC(x_k,x_new,xbest,X)
    eval=1;
    LB=min(X); UB=max(X);
    dbk=di(x_k,xbest,LB,UB);
    db1=di(x_new,xbest,LB,UB);
 
    if db1 >= dbk,
        eval=0;
    end
end