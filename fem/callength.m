% Calculate element length
function len = callength(Node,Ele,NoE,nvar)
    global L H
    H = zeros(nvar,NoE); L = zeros(1,NoE);
    for e = 1:NoE
        C = [Node(Ele(e,1),1), Node(Ele(e,1),2), Node(Ele(e,1),3); 
             Node(Ele(e,2),1), Node(Ele(e,2),2), Node(Ele(e,2),3)];
        x1 = C(1,1); y1 = C(1,2); z1 = C(1,3);
        x2 = C(2,1); y2 = C(2,2); z2 = C(2,3);
        dx = (x2-x1); dy = (y2-y1); dz = (z2-z1);        
        L(e) = sqrt(dx^2 + dy^2 + dz^2);
        
        H(Ele(e,3),e) = 1;
    end    
    len = L;
end