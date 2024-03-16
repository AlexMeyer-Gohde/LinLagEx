%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%input_example_me_3.m
%
%This program:
%               is the model file used for estimation with the truncated 
%               model. I.e., used in calculating the log-likelihood function 
%               for a given set of parameter found in the vector parameters.
%
%           Note the values under the heading options:we don't want to
%           calculate anything else other than the recursive solution. 
%
%           Note also the mapping from the vector "parameters" to the model
%           parameters.
%
%THIS VERSION: 1.1 December 9, 2009
%
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Paramters
lambda=parameters(1);
xi=parameters(2);
a_1=1/parameters(3);
phi_pi=parameters(4);
phi_y=parameters(5);
rho_is=parameters(6);
rho_pc=parameters(7);
phi_R=parameters(8);
sigma_eps_is=parameters(9);
sigma_eps_pc=parameters(10);
sigma_eps_mp=parameters(11);





%Variable names
ENDOGENOUS_VARIABLE_NAMES={'Inflation'
'Output Gap'
'Nominal Interest Rate'
'PC'};
EXOGENOUS_VARIABLE_NAMES={'IS'
'PC'
'MP'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0=sum A_iE_t-i(x_t+1)+sum B_iE_t-i(x_t)+sum C_iE_t-i(x_t-1)+sum
% F_iE_t-i(z_t+1)+sumG_iE_t-i(z_t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matrices

A_0=[0 0 0  0
    a_1 1 0  0
    0 0 0  0
    0 0 0  0];
B_0=[-1 (1-lambda)*xi/lambda 0  0
    0 -1 -a_1  0
(1-phi_R)*phi_pi   (1-phi_R)*phi_y   -1  0
0 0 0  -1];
C_0=[zeros(2,4) 
        0 0 phi_R  0
        0 0 0  0];
F_0=[zeros(4,3)];
G_0=[0 (1-lambda)/lambda 0
    1 0 0
0 0 1
0 1 0];




%counter=maximum_counter_value
%name of the counter used in subsequent matrices
it_name='J_1';
%maximum value of counter: either 'infinity' or a scalar (without
%quotation marks)
it_max_value=3;


%counter matrices

A_j='[zeros(4,4)]';
B_j='[(1-lambda)*lambda^(J_1-1) xi*(1-lambda)*lambda^(J_1-1) 0  (1-lambda)*lambda^(J_1-1); zeros(3,4)]';
C_j='[0 -xi*(1-lambda)*lambda^(J_1-1) 0  -(1-lambda)*lambda^(J_1-1); zeros(3,4)]';
F_j='[zeros(4,3)]';
G_j='[zeros(4,3)]';

%limit matrices
%If a counter goes to infinity, please enter the limiting summed matrices
clear A_inf B_inf C_inf F_inf G_inf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define exogenous process as a VAR(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AR-coefficients
N=[rho_is,0,0;0,rho_pc,0;0,0,0];
%Covariance matrix
Omega=[sigma_eps_is^2 0 0
0 sigma_eps_pc^2 0
0 0 sigma_eps_mp^2];

%Options
PERIOD=4;
tolerance=1e-5;
RUN_IMPULSE=0;
RUN_ANTICIPATED_IMPULSE=0;
RUN_SIMULATION=0;
RUN_SPECTRAL=0;
solution_method='AIM';
linlagex;