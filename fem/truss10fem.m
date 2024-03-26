% ----------------------------------------------------------------------- %
% FEM of the 10-bar truss optimization 
% Pham Hoang Anh, 2024
% Department of Structural Mechanics, 
% Hanoi University of Civil Engineering
% Email: anhph2@huce.edu.com
% ------------------------------------------------------------------------%
function [D, S] = truss10fem(x)

global E 
global NoN NoE Node Ele Nload fixdofs
NoD = NoN*2;

%% Element conectivity   
EleID = zeros(NoE,4);
EleID(1:NoE,:) = [2*Ele(1:NoE,1)-1, 2*Ele(1:NoE,1), 2*Ele(1:NoE,2)-1, 2*Ele(1:NoE,2)];
iK = reshape(kron(EleID,ones(4,1))',16*NoE,1);
jK = reshape(kron(EleID,ones(1,4))',16*NoE,1);

%% Caculate element stiffness matrix 
sK = zeros(16,NoE); 
Tes = zeros(2,4,NoE);
Les = zeros(NoE,1);
for i=1:NoE,
    xi = Node(Ele(i,1:2),1); yi = Node(Ele(i,1:2),2);
    Li = sqrt(((xi(2)-xi(1))^2 + (yi(2)-yi(1))^2));
    Ai = x(Ele(i,3));
    % stiffness matrix
    Ki = E*Ai/Li*[1 -1; -1  1];
    
    % transform matrix   
    Co = (xi(2)-xi(1))/Li; Si = (yi(2)-yi(1))/Li;
    Te = [Co Si   0      0 ;
          0    0   Co   Si];    
    Tes(:,:,i) = Te;
    Les(i) = Li;
    Ke = (Te')*Ki*Te; sK(:,i) = Ke(:);
end

%% Assemble for global stiffness matrix and load vector
sK=sK(:);
K = sparse(iK,jK,sK); 

% Add loading
F = zeros(NoD,1);
for i=1:size(Nload,1)
    node = Nload(i,1);
    F(node*2-1:node*2) = Nload(i,2:3);
end

%% Calculate displacement
alldofs = 1:NoD;
D = zeros(NoD,1);
freedofs = setdiff(alldofs,fixdofs);
D(freedofs,1) = K(freedofs,freedofs)\F(freedofs);

%% Calculate stress
S = zeros(NoE,1);
for i=1:NoE
    De = Tes(:,:,i)*D(EleID(i,:));
    dL = De(2) - De(1);
    S(i) = E*dL/Les(i);
end

end
