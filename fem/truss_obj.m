% Calculate truss weight
function f = truss_obj(x)
    global ro L H
    A = x*H;   
    f = sum(ro*(L.*A));
end