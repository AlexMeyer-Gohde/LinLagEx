%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%input_example_meyer_gohde_2009.m
%
%This is the model
%from: 
%Meyer-Gohde, A. (2009): "Linear Rational Expectations Models with Lagged
%Expectations: A Synthetic Method
%
%THIS VERSION: 1.1 December 9, 2009
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('..\linlagex')

%Paramters


lambda=0.6875;
xi=0.2264;
a_1=1/3.3995;
phi_pi=1.3238;
phi_y=0.2487;
rho_is=0.8839;
rho_pc=0.7649;
phi_R=0.7641;
sigma_eps_is=0.1614;
sigma_eps_pc=0.6609;
sigma_eps_mp=0.2656;





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
it_max_value='infinity';


%counter matrices

A_j='[zeros(4,4)]';
B_j='[(1-lambda)*lambda^(J_1-1) xi*(1-lambda)*lambda^(J_1-1) 0  (1-lambda)*lambda^(J_1-1); zeros(3,4)]';
C_j='[0 -xi*(1-lambda)*lambda^(J_1-1) 0  -(1-lambda)*lambda^(J_1-1); zeros(3,4)]';
F_j='[zeros(4,3)]';
G_j='[zeros(4,3)]';

%limit matrices
%If a counter goes to infinity, please enter the limiting summed matrices
A_inf=[A_0];
B_inf=[0 xi/lambda 0  1
    0 -1 -a_1  0
(1-phi_R)*phi_pi   (1-phi_R)*phi_y   -1  0
  0 0 0  -1];
C_inf=[0 -xi 0  -1
        0 0 0  0 
        0 0 phi_R 0
        0 0 0  0];
F_inf=[F_0];
G_inf=[G_0];

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
linlagex;