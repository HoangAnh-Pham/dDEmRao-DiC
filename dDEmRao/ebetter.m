%% e constrained comparison
function ecom = ebetter(f1,con1,f2,con2,e)
    ecom = 0;
    if ((con1 <= e)&&(con2 <= e)) || (con1 == con2), 
        if f1 < f2, ecom = 1; end;
    else if con1 < con2, ecom = 1; end;
    end
end