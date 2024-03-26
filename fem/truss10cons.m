% Calculate constraint violation
function [C,Ceq] = truss10cons(x)

    s_max = 25;     % Stress limit  (ksi)
    d_max = 2;      % Displacement limit  (in)

    % FEM
    [D,S] = truss10fem(x);
    sc = abs(S); 
    dc = abs(D);

    C1 = sc/s_max - 1; 
    C2 = dc/d_max - 1; 
    C = [C1;C2];
    Ceq = [];
end