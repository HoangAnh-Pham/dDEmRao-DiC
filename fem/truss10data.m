% ----------------------------------------------------------------------- %
% FEM data of the 10-bar truss optimization 
% Hoang-Anh Pham, 2024
% Department of Structural Mechanics, 
% Hanoi University of Civil Engineering
% Email: anhph2@huce.edu.com
% ------------------------------------------------------------------------%
function truss10data()

global nvars XB 
global E ro
global NoN NoE Node Ele Nload fixdofs

nvars = 10;     % Number of design variables
XB = [0.1,40];  % Variable bounds

E = 1.0e+4;     % Elastic modulus [ksi]
ro = 0.1e-0;    % Material density [lb/in3]

%% Node data: [x-coord, y-coord, z-coord]
% Node coordinate
a = 360.;               % Span length [in]
% Node(i) = [xi, yi, zi]
Node = [2*a,a,0; 
        2*a,0,0;
        a, a, 0; 
        a, 0, 0
        0, a, 0; 
        0, 0, 0];
NoN = size(Node,1);     % number of nodes

% Restraints
fixdofs = [9 10 11 12];

%% Element data: [node1, node2, sectionID]
Ele = [5, 3, 1
       3, 1, 2
       6, 4, 3
       4, 2, 4
       3, 4, 5
       1, 2, 6
       5, 4, 7
       6, 3, 8
       3, 2, 9
       4, 1, 10];
NoE = size(Ele,1);      % number of elemenents

%% Load condition
P1 = -100.;      % Concentrated load [kips]
P2 = 0.;
% Nodal load: [node, Fx, Fy, Fz]
Nload = [1 0 P2; 
         2 0 P1; 
         3 0 P2; 
         4 0 P1];
     
%% Calculate element length     
callength(Node,Ele,NoE,nvars);

end
