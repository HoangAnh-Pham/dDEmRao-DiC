% ------------------------------------------------------------------------%
% MAIN PROGRAM
% Enhanced Differential Evolution-Rao with Distance Comparison Method
% <dDEmRao>
% Hoang-Anh Pham, 2024
% Department of Structural Mechanics, 
% Hanoi University of Civil Engineering
% Email: anhph2@huce.edu.com
% ------------------------------------------------------------------------%
addpath('fem');
addpath('dDEmRao');

clc; close all; clear all; 
global nvars XB 
%% Setting optimization parameters
Ng = 1000;      % No. iterations
Tol = -1e-6;    % Min. relative error
NP = 25;        % Population size
NoR = 1;        % No. runs
para = [Ng, Tol, NP];

%% Optimization problem
problem = {{'3-bar',@truss3obj,@truss3cons,@truss3data};        % 1
           {'10-bar',@truss_obj,@truss10cons,@truss10data}};    % 2

pb = 2;  % problem ID
truss_name = problem{pb}{1};
fname = problem{pb}{2};
fcons = problem{pb}{3};
data = problem{pb}{4};

feval(data);

LB = min(XB)*ones(1,nvars); % Lower bound
UB = max(XB)*ones(1,nvars); % Upper bound
DX = [];

disp(['Problem: ',truss_name]);

%% Optimization method
method={{'DE',          @dDEmRao,'rnd','',''};          % 1        
        {'Rao1',        @dDEmRao,'rao','',''};          % 2
        {'DE-Rao1',     @dDEmRao,'hb1','',''};          % 3
        {'DE-Rao1-DiC', @dDEmRao,'hb1','','dic'};       % 4
        {'dDEmRao-DiC', @dDEmRao,'hb2','d','dic'}};     % 5
        
ID = [1,5]; % List of method ID
method = method(ID); 
mt = 1:length(method);

%% Run optimization
for i=1:length(mt)
    tic;
    disp(['Method: ',method{mt(i)}{1}]);
    algorithm = method{mt(i)}{2};
    option = {method{mt(i)}{3:end}};
    varin = {algorithm,fname,fcons,nvars,LB,UB,DX,para,option{:}};
    
    for t=1:NoR
        disp(['Run: ',int2str(t)]);
        [xopt,fopt,exitflag,out,X,scores,V,FE,DI,S] = feval(varin{:});
        
        % show results
        disp(['Optimized weight:',num2str(fopt)]);
        disp(['CV:',num2str(max(feval(fcons,xopt)))]);
        disp(out);
        
        figure; hold all; box on;
        title([truss_name,', ',method{mt(i)}{1},' ','run ',num2str(t)]);
        plot(FE(2,:),'-b','LineWidth',1.5); 
        plot(FE(3,:),'-r','LineWidth',1.5); 
        
        ylabel('No. generated solutions');
        xlabel('Iterations');
        legend('by Rao','by DE');
        hold off;
        
        figure; box on;
        plot(V,'LineWidth',1.5); 
        title([truss_name,', ',method{mt(i)}{1},' ','run ',num2str(t)]);
        xlabel('FEs'); ylabel('Weight');
        hold off;
        
        figure; box on;
        semilogy(DI,'LineWidth',1.5); 
        title([truss_name,', ',method{mt(i)}{1},' ','run ',num2str(t)]);
        xlabel('Iterations'); ylabel('Diversity index'); box on;
        hold off;
        
        % save results
        save([truss_name,'-',method{mt(i)}{1},'_T',num2str(t),'.mat']);
    end
    toc;
end


