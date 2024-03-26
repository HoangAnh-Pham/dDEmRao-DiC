% ------------------------------------------------------------------------%
% Enhanced Differential Evolution-Rao with Distance Comparison Method
% <dDEmRao>
% Hoang-Anh Pham, 2024
% Department of Structural Mechanics, 
% Hanoi University of Civil Engineering
% Email: anhph2@huce.edu.com
% ------------------------------------------------------------------------%
function [xopt,fopt,exitflag,output,X,Fitness,history,FE,DI,S] = dDEmRao(fname,fcons,nvars,LB,UB,DX,option,type,dir,DiC)

    % Setting algorithm parameters
    T_max = option(1); 
    Tol = option(2);
    NP = option(3);
    Fk = 0.8; 
    CR = 0.9; 

    % Initialize the population/solutions
    [row,col] = size(DX);
    rnd = rand(NP,nvars);
    lb = repmat(LB,NP,1);
    ub = repmat(UB,NP,1);
    X = lb + rnd.*(ub-lb);
    for i=1:NP
        for j=1:row, X(i,DX(j,1)) = dismap(X(i,DX(j,1)),DX(j,2:col)); end     
        % Evaluate fitness and constraint violation
        Fitness(i) = feval(fname,X(i,:)); 
        Const(i) = constraint(fcons,X(i,:),1e-4); 
    end    
    
    % Sort the population in ascending order of fitness       
    [~,ID] = consort(Fitness,Const,20,0); 
    xbest = X(ID(1),:); 
    fbest = Fitness(ID(1)); 
    cbest = Const(ID(1));    
    fmean = mean(Fitness);
    e0=Const(ID(round(0.2*NP))); %e0=0;
    
    % Calculate diversity index of initial population
    DI_0 = d_index(X,LB,UB); 
    DI_t = DI_0;
    DI = DI_t; 
    
    % Initializing 
    num_FE = NP; num_succ = 0; num_fail = 0; 
    num_skip = 0; 
    S = num_skip; 
    
    N = [0,0]; history = []; 
    FE = num_FE;
    no_DE=0; no_Rao=0;
    FE_DE=no_DE; FE_Rao=no_Rao;
    
%    fprintf('\n						 Best\t\t    Mean\n')
%    fprintf('Generation\tNFE\t\t f(x)\n')
    
    % Start the iterations
    iter=0; dev = Tol+1;
    while (iter < T_max) && (dev > Tol),
       
%       fprintf('\t%i\t\t%i\t\t\t %f\t %f\t\n',...
%            iter,num_FE,fbest,fmean)

       iter=iter+1;
       
       inum_FE=0; inum_skip = 0; 
       ino_DE=0; ino_Rao=0;

       en = e0*(1+(DI_t-DI_0)/DI_0)^2; if en<0, en=0; end

       % Select method
       switch type
          % Rao-1
          case 'rao', DR=2;
          % DE/rand/1
          case 'rnd', DR=-1;
          % DE-Rao1
          case 'hb1', DR = DI_t/DI_0; 
          % dDEmRao
          case 'hb2', DR = DI_t/DI_0; 
       end          
        
       % Sort the population in ascending order of fitness       
       [~,ID] = consort(Fitness,Const,20,0); 
       
       xbest = X(ID(1),:);
       fbest = Fitness(ID(1)); 
       cbest = Const(ID(1)); 
       xold = xbest; 
       fold = fbest; 
       xworst = X(ID(end),:);

       X = X(ID(1:end),:); 
       Fitness = Fitness(ID(1:end));
       Const = Const(ID(1:end));
       X1=X; Fitness1=Fitness; Const1=Const; 
              
       % Loop over all solutions
       for k=1:NP
           % Take base vector from current solutions
           JK = randperm(round(NP)); JK(ID(JK)==k)=[]; 
           r1 = ID(JK(1)); 

           % Select randomly two other vectors for mutation
           r2 = k; while (r2==k)||(r2==r1), r2 = randperm(NP,1); end;
           r3 = k; while (r3==k)||(r3==r1)||(r3==r2), r3 = randperm(NP,1); end;

           % Set direction for mutation
           if dir=='d'
                ecom = ebetter(Fitness(r2),Const(r2),Fitness(r3),Const(r3),en); 
                if ecom==1, d=1; else d=-1; end;   
           else d=1;
           end
           
           %% Hybrid Rao - DE
           if rand^2>DR,
                % DE mutation and crossover;        
                x_new = X(k,:);
                xbase = X(r1,:);
                r = randperm(nvars,1);
                for i = 1:nvars               
                    if (rand <= CR) || (i == r), 
                        x_new(i) = xbase(i) + d*Fk*(X(r2,i)-X(r3,i));       
                    end
                end
                ino_DE=ino_DE+1;
           else
                % Rao-1
                x_new = X(k,:) + rand(1,nvars).*(xold-xworst);
                if strcmp(type,'hb2' ),
                    % Modified Rao-1 operator
                    x_new = x_new + d*rand(1,nvars).*(X(r2,:)-X(r3,:));
                end
                ino_Rao=ino_Rao+1;
           end
  
           %%
           % Check bound constraints
           x_new = boundConst(x_new,X(k,:),LB,UB);
           for j=1:row, x_new(DX(j,1)) = dismap(x_new(DX(j,1)),DX(j,2:col)); end     
           
           % Apply the Distance Comparison
           eval=1; 
           switch DiC  
                case 'dic' 
                    eval=DIC(X(k,:),x_new,xbest,X); 
                    %if Const(k)>en || k==1, eval=1; end
                    if Const(k)>en, eval=1; end                 
           end
           if eval==0,
              inum_skip=inum_skip+1;
           end
           
           % Check to skip usefuless evaluation
           if eval>0, eval=0;
              fnew = feval(fname,x_new); 
              
              inum_FE = inum_FE + 1; 

              if Const(k)>en, eval=1;
              else if fnew<Fitness(k), eval=1; end
              end
              history = [history,fbest];
           end
         
           if eval>0,           
                cnew = constraint(fcons,x_new,1e-4);
                
                % Select new solution as member if better
                ecom  = ebetter(fnew,cnew,Fitness(k),Const(k),en); 
                if (ecom == 1), 
                    X1(k,:) = x_new; 
                    Fitness1(k) = fnew; 
                    Const1(k) = cnew;
                    X(k,:) = x_new; 
                    Fitness(k) = fnew; 
                    Const(k) = cnew;
                    
                    num_succ = num_succ + 1; 
                else
                    num_fail = num_fail + 1;
                end;                
           end      
           
           % Update the current global best
           ecom = ebetter(Fitness1(k),Const1(k),fbest,cbest,0);
           if ecom==1, xbest = X1(k,:); fbest = Fitness1(k); cbest = Const1(k); end                          
           if fbest<fold, N=[iter,num_FE]; fold=fbest; end
           
       end % End for
       X=X1; Fitness=Fitness1; Const=Const1;
       
       DI_t = d_index(X,LB,UB); DI=[DI,DI_t];
       no_DE=no_DE + ino_DE;
       no_Rao=no_Rao + ino_Rao;
       num_FE = num_FE + inum_FE; 
       num_skip = num_skip + inum_skip;
       
       S = [S,num_skip/(iter*NP)];
       FE = [FE,num_FE];
       FE_DE=[FE_DE,no_DE]; 
       FE_Rao=[FE_Rao,no_Rao];
       
       fmean = mean(Fitness); 
       dev = abs(fmean/fbest-1);
      
    end % End while
    
    % Return the optimal solution
    xopt = xbest; fopt = fbest;
    FE = [FE; FE_Rao; FE_DE];
    
    if iter>=T_max, 
        exitflag = 0;   massage = 'generations exceed max. generation';
    else exitflag = 1;  massage = 'diversity is lower than Tol';
    end
    
    output = struct('generations',iter,...
                    'skiprate',num_skip/(num_FE+num_skip)*100,...
                    'succskip',0,...                    
                    'funccount',num_FE,...
                    'constcount',num_succ+num_fail,...
                    'success',num_succ,...
                    'fail',num_fail,...
                    'minfunc',N,...
                    'maxconstraint',cbest,...
                    'diversity',d_index(X,LB,UB),...
                    'message',massage);    
    
end % Main function



