%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%input_example_wang_wen_1_2006.m
%
%This is the first example 
%(Sticky Information)
%from: 
%Wang, P., and Y. Wen (2006): “Solving Linear Difference Systems with
%Lagged Expectations by a Method of Undetermined Coefficients,” Working
%Papers 2006-003, Federal Reserve Bank of St. Louis.
%
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
alpha=0.2;  % capital share
ela_p=10;   % substitution elasticity of intermediate goods
delta=0.025;% capital depreciate rate
beta=0.99;  % discount rate
theta=0.8;  % degree of nominal rigidity
gamac=0.005;% risk aversion parameter in consumption 
gaman=0.05; % inverse of labor supply elasticity
rho=0.6;    % shock persistence %
 
mc=(ela_p-1)/ela_p;    % steady-state marginal cost
yk=(1+(delta-1)*beta)/(beta*mc)/alpha; % steady-state yk ratio
ck=yk-delta;

%Variable names
ENDOGENOUS_VARIABLE_NAMES={'Consumption'
'Inflation'
'Labor'
'Capital'
'Marginal Costs'
'Output'
'Investment'
'Nominal Interest Rate'};
EXOGENOUS_VARIABLE_NAMES={'Money Growth'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0=sum A_iE_t-i(x_t+1)+sum B_iE_t-i(x_t)+sum C_iE_t-i(x_t-1)+sum
% F_iE_t-i(z_t+1)+sumG_iE_t-i(z_t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matrices
%       c   pi  n   k   mc  y   i   R
A_0=[   gamac 1 0   0   0   0   0   0
        -gamac*beta*(1-delta) 0 (1-beta*(1-delta))*(gaman+1) 0 0 0 0 0
         0  0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   1   -(gaman+alpha)  0   1   0   0   0];


B_0=[   0   0  (gaman+alpha)    0   -1  0   0   0
        gamac   0   0  -(1-beta*(1-delta))  0   0   0   0 
        0   -1  0   0   (1-theta)/theta 0   0   0
        ck  0   0   0   0   -yk  delta 0
        0   0   0   -1  0   0   delta   0
        0   0   (1-alpha)   0   0   -1  0   0
        0   -1  0   0   0   -1  0   0
        0   0   (gaman+alpha)   alpha   -1  0   0   -1];

C_0=[   0   0   0   -alpha  0   0   0   0
        zeros(3,8)
        0   0   0   (1-delta)   0   0   0   0
        0   0   0   alpha   0   0   0   0   
        0   0   0   0   0   1   0   0   
        0   0   0   -alpha  0   0   0   0];
F_0=[zeros(8,1)];
G_0=[zeros(6,1)
    1
    0];


%name of the counter used in subsequent matrices
it_name='J_1';
%maximum value of counter: either 'infinity' or a scalar (without
%quotation marks)
it_max_value='infinity';


%counter matrices
%       c   pi  n   k   mc  y   i   R
A_j='[zeros(8,8)]';
B_j='[zeros(2,8); 0 (1-theta)*theta^(J_1-1) 0 0 (1-theta)*theta^(J_1-1) 0 0 0; zeros(5,8)]';
C_j='[zeros(2,8); 0 0 0 0 -(1-theta)*theta^(J_1-1) 0 0 0; zeros(5,8)]';
F_j='[zeros(8,1)]';
G_j='[zeros(8,1)]';

%limit matrices
%If a counter goes to infinity, please enter the limiting summed matrices
%Alternatively, the limits can be calculated symbolically. 
A_inf=[A_0];
B_inf=B_0+[zeros(2,8); 0 1 0 0 1 0 0 0; zeros(5,8)];
C_inf=C_0+[zeros(2,8); 0 0 0 0 -1 0 0 0; zeros(5,8)];
F_inf=[F_0];
G_inf=[G_0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define exogenous process as a VAR(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AR-coefficients
N=[rho];
%Covariance matrix
Omega=[1];

%Options
PERIOD=4;


linlagex;
