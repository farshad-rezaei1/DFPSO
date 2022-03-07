%___________________________________________________________________%
% Dual Fitness Particle Swarm Optimization (DFPSO) Algorithm        %
%                                                                   %
% Developed in MATLAB R2018b                                        %
%                                                                   %
% Inventor and programmer: Farshad Rezaei, PhD                      %
%                                                                   %
% e-Mail: farshad.rezaei@gmail.com                                  %
%         f.rezaei@alumni.iut.ac.ir                                 %
%                                                                   %
% Homepage: https://www.linkedin.com/in/farshad-rezaei-5a92559a/    %
%                                                                   %
% Main paper: Rezaei, F., Safavi, H.R. Sustainable Conjunctive      %
% Water Use Modeling Using Dual Fitness Particle Swarm Optimization %
% Algorithm. Water Resour Manage 36, 989–1006 (2022).               %
% https://doi.org/10.1007/s11269-022-03064-w                        %
%___________________________________________________________________%

% DFPSO algorithm                                                                  
function [z_iter,z_final,pos_final] = DFPSO(np,nx,maxit,varmax,varmin,velmax,velmin,epsilon,k_max,k_min,fobj)
% disp(['Number of Iterations = ',num2str(it)]);
pp_pbest=zeros(np,nx);
DFI=zeros(np);
optimal_pos=zeros(1,nx);
z_pbest=zeros(np);
z_optimal=inf*ones(maxit);
pos_final=zeros(nx);
z_iter=zeros(maxit);

% Initialization process of the algorithm
[pp,pv]=Initialization(np,nx,varmax,varmin,velmax,velmin);

% Start the optimization process
it=1;

% Objective function evaluations and determine the personal best solutions and objectives
for j=1:np
    z=fobj(pp(j,1:nx));
    z_pbest(j)=z;
    pp_pbest(j,1:nx)=pp(j,1:nx);
end

% Calculating Mean and Std of the Pbests' objective values
ave=mean(z_pbest(1:np));
stdev=std(z_pbest(1:np));

% Evaluating the Dual Fitness Index for the Particles
prod=ones(np);
for j=1:np
    for jj=1:np
        if jj~=j
            prod(j)=prod(j)*(1/(1+exp((-4)/(stdev*sqrt(exp(1)))*(z_pbest(jj)-ave)))); % Eq.(8)
        end
    end
end
sum_prod=sum(prod(1:np));

% Calculating the unique Global best guide for each particle
pp_gbest=zeros(np,nx);
for j=1:np 
    for jj=1:np
        if jj~=j
            DFI(jj)=prod(jj)/(sum_prod+epsilon);
            pp_gbest(j,1:nx)=pp_gbest(j,1:nx)+DFI(jj)*pp_pbest(jj,1:nx); % Eq.(9)
        end
    end
end

% Determining the best-so-far objective value and the best-so-far particle
for j=1:np
    if z_pbest(j)<z_optimal(it)
        z_optimal(it)=z_pbest(j);
        optimal_pos(it,:)=pp_pbest(j,:);
    end
end

% Save the best-so-far objective value in the current run
z_iter(it)=z_optimal(it);

% The Main Loop
while it<maxit
    it=it+1;
    k=k_max-(k_max-k_min)*(it/maxit); % Eq.(5)
%     disp(['Number of Iterations = ',num2str(it)]);
    for j=1:np 
        phi1=2*rand(1,nx);phi2=2*rand(1,nx);phi=phi1+phi2;
        khi=(2*k)./abs(2-phi-sqrt(phi.*(phi-4))); % Eq.(4)
        
        % Update the velocity of the particles- Eq.(10)
        pv(j,1:nx)=khi.*(pv(j,1:nx)+phi1.*(pp_pbest(j,1:nx)-pp(j,1:nx))+...
            phi2.*(pp_gbest(j,1:nx)-pp(j,1:nx)));
                
        % Return back the velocity of the particles if going beyond the velocity boundaries
        flag4lbv=pv(j,:)<velmin(1,:);
        flag4ubv=pv(j,:)>velmax(1,:);
        pv(j,:)=(pv(j,:)).*(~(flag4lbv+flag4ubv))+velmin.*flag4lbv+velmax.*flag4ubv;
        
        % Update the position of the particles- Eq.(2)
        pp(j,:)=pp(j,:)+pv(j,:);
        
        % Return back the position and velocity of the particles if going beyond the position boundaries
        flag4lbp=pp(j,:)<varmin(1,:);
        flag4ubp=pp(j,:)>varmax(1,:);
        pp(j,:)=(pp(j,:)).*(~(flag4lbp+flag4ubp))+varmin.*flag4lbp+varmax.*flag4ubp; 
        pv(j,:)=(pv(j,:)).*(ones(1,nx)-2*(flag4lbp+flag4ubp)); 
    
        % Objective function evaluations and determining the personal best solutions and objectives
        z=fobj(pp(j,1:nx));
        if z<z_pbest(j)
            z_pbest(j)=z;
            pp_pbest(j,1:nx)=pp(j,1:nx);
        end
    end
    
    % Calculating Mean and Std of the Pbests' objective values
    ave=mean(z_pbest(1:np));
    stdev=std(z_pbest(1:np));

    % Evaluating the Dual Fitness Index for the Particles
    prod=ones(np);
    for j=1:np
        for jj=1:np
            if jj~=j
                prod(j)=prod(j)*(1/(1+exp((-4)/(stdev*sqrt(exp(1)))*(z_pbest(jj)-ave)))); % Eq.(8)
            end
        end
    end
    sum_prod=sum(prod(1:np));
    
    % Calculating the unique Global best guide for each particle
    pp_gbest=zeros(np,nx);
    for j=1:np 
        for jj=1:np
            if jj~=j
                DFI(jj)=prod(jj)/(sum_prod+epsilon);
                pp_gbest(j,1:nx)=pp_gbest(j,1:nx)+DFI(jj)*pp_pbest(jj,1:nx); % Eq.(9)
            end
        end
    end
    
    % Determining the best-so-far objective value and the best-so-far particle
    for j=1:np
        if z_pbest(j)<z_optimal(it)
            z_optimal(it)=z_pbest(j);
            optimal_pos(it,:)=pp_pbest(j,:);
        end
    end
    
    % Save the best-so-far objective value in the current run
    z_iter(it)=z_optimal(it);
end

% Save the final best solution and objective revealed upon the end of the optimization process
z_final=z_optimal(maxit);
pos_final(1:nx)=optimal_pos(maxit,1:nx);
end