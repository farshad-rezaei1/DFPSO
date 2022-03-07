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

% The initial parameters that you need are:
%__________________________________________
% fobj = @YourCostFunction
% nx = number of your variables
% lb= the lower bound of variables which can generally be a fixed number or a vector
% ub= the upper bound of variables which can generally be a fixed number or a vector
% notice: if the lower and upper bounds are not fixed for all variables, 
% they appear in the forms of the vectors "varmin" and "varmax", as illustrated in following

% To run DFPSO: [z_iter,z_final,pos_final]=DFPSO(np,nx,maxit,varmax,varmin,velmax,velmin,epsilon,k_max,k_min,fobj);
%__________________________________________
% Set the required parameters to run the DFPSO algorithm

% This code is for solving the minimization problems. To maximize a desired 
% cost function,please implement this code upon inverting the cost function

clc
clear
close all
tic
run=1; % Maximum number of the algorithm runnings conducted
np=50; % Number of search agents
maxit=1000; % Maximum number of iterations
Function_name='F1'; % Name of the test function that can be from F1 to F23
epsilon=eps; % This parameter must be either set to eps (typically for the uni-modal functions) or zero (typically for the multi-modal or complex functions)
k_max=0.9; % Upper bound of the inertia weight
k_min=0.4; % Lower bound of the inertia weight
[lb,ub,nx,fobj]=Objective_Function(Function_name); % Load details of the selected benchmark function
varmax=ub*ones(1,nx); % Upper bound defined for the positions which can generally be a desired vector
varmin=lb*ones(1,nx); % Lower bound defined for the positions which can generally be a desired vector
limvel=0.1; % A ratio of the maximum distance in the search space to form the maximum velocity 
velmax=limvel*(varmax(1,1:nx)-varmin(1,1:nx)); % Upper bound defined for the velocities
velmin=-velmax; % Lower bound defined for the velocities
z_iter_main=zeros(run,maxit);
z_final_main=zeros(run);
pos_final_main=zeros(run,nx);
x1=zeros(maxit);
y1=zeros(maxit);

% Run the DFPSO algorithm for "run" times 
for nrun=1:run
    [z_iter,z_final,pos_final]=DFPSO(np,nx,maxit,varmax,varmin,velmax,velmin,epsilon,k_max,k_min,fobj);
     z_iter_main(nrun,1:maxit)=z_iter(1:maxit);
     z_final_main(nrun)=z_final;
     pos_final_main(nrun,1:nx)=pos_final(1:nx);
%      disp(['The best objective function value obtained by DFPSO = ',num2str(z_final_main(nrun))]);
%      disp(['The best solution obtained by DFPSO = ','[',num2str(pos_final_main(nrun,1:nx)),']']);
end

% Display the comprehensive results
disp(['The final statistical results calculated when implementing the DFPSO algorithm for ',num2str(run),' times are as follows:']);
disp(['The average of the final objective function values calculated over ',num2str(run),' times = ',num2str(mean(z_final_main(1:run)))]);
disp(['The median of the final objective function values calculated over ',num2str(run),' times = ',num2str(median(z_final_main(1:run)))]);
disp(['The best of the final objective function values calculated over ',num2str(run),' times = ',num2str(min(z_final_main(1:run)))]);
disp(['The standard deviation of the final objective function values calculated over ',num2str(run),' times = ',num2str(std(z_final_main(1:run)))]);

% Plot the chart of the objective function values obtained by DFPSO over the course of iterations
for i=1:maxit
    x1(i)=i;sum1=0;
    for j=1:run
        sum1=sum1+z_iter_main(j,i);
    end
    y1(i)=sum1/run;
end
semilogy(x1,y1,'-r')
xlabel('Iteration');
ylabel('Average best-so-far');
legend('DFPSO');
hold on
time_dfpso = toc;
disp(['Elapsed time of running the DFPSO for ',num2str(run),' times = ',num2str(time_dfpso),' seconds']);