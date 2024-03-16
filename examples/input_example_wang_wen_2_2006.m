%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%input_example_wang_wen_2_2006.m
%
%This is the second example 
%(Labor Hoarding)
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
T=1369;  
N=50;
xi=60;   
f=324.77754855;
beta=0.99;  
alpha=0.36;  
rho_theta=0.9;
rho_g=0.9;
rho_a=0.9;
u_bar=1;

s_g=0.2;
delta_bar=0.025;


e_bar=fzero(@(x)log(T)-log(T-xi-x*f)-f*x/(T-xi-x*f),1);
pi=e_bar*f/(T-xi-e_bar*f);
eta=(1-beta*(1-delta_bar))*(1-alpha);
phi=(1-beta*(1-delta_bar))/(beta*delta_bar);
delta=delta_bar/(u_bar^phi);
k_y_bar=beta*alpha/(1-beta*(1-delta_bar));
s_i=delta_bar*(k_y_bar);
s_c=1-s_i-s_g;



%Variable names
ENDOGENOUS_VARIABLE_NAMES={'Costate'
'Utilization'
'Capital'
'Effort'
'Employment'
'Consumption'
'Output'};
EXOGENOUS_VARIABLE_NAMES={'Productivity'
'Consumption Urge'
'Government Spending'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0=sum A_iE_t-i(x_t+1)+sum B_iE_t-i(x_t)+sum C_iE_t-i(x_t-1)+sum
% F_iE_t-i(z_t+1)+sumG_iE_t-i(z_t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matrices
%       L   u   k   e   n   c   y
A_0=[   zeros(4,7);
        1 -eta 0 eta eta 0   0
        zeros(2,7)];

%       L   u   k   e   n   c   y
B_0=[   0   0   0   0   -alpha 0 0
        1   0   0   0   0   1   0
        0   -(phi-alpha) 0 (1-alpha) (1-alpha) 0 0
        1   alpha 0 -(pi+alpha) -alpha 0 0
        -1 0    -eta   0   0   0   0
        0   (alpha-s_i*phi) -s_i/delta_bar (1-alpha) (1-alpha) -s_c 0
        0   alpha 0 (1-alpha) (1-alpha) 0 -1];

C_0=[   zeros(2,7)
        0   0   -(1-alpha)  0   0   0   0
        0   0   alpha       0   0   0   0
        zeros(1,7)
        0   0   (alpha+s_i*(1-delta_bar)/delta_bar) 0 0 0 0
        0   0   alpha   0   0   0   0];
%       A   th  g
F_0=[zeros(4,3)
        eta/(1-alpha) 0 0
    zeros(2,3)];
G_0=[   0   0   0
        0   -1  0
        1   0   0
        1   0   0
        0   0   0
        1   0   -s_g
        1   0   0];


%name of the counter used in subsequent matrices
it_name='J_1';
%maximum value of counter: either 'infinity' or a scalar (without
%quotation marks)
it_max_value=N;


%counter matrices
%       L   u   k   e   n   c   y
A_j='[zeros(7,7)]';
B_j='[zeros(7,7)]+(J_1==it_max_value)*[1 alpha 0 -alpha 0 0 0;zeros(6,7)]';
C_j='[zeros(7,7)]+(J_1==it_max_value)*[0 0 alpha 0 0 0 0;zeros(6,7)]';
F_j='[zeros(7,3)]';
G_j='[zeros(7,3)]+(J_1==it_max_value)*[1  0 0;zeros(6,3)]';

%limit matrices
%If a counter goes to infinity, please enter the limiting summed matrices
%Alternatively, the limits can be calculated symbolically. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define exogenous process as a VAR(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AR-coefficients
N=[rho_a 0 0
    0 rho_theta 0
    0 0 rho_g];
%Covariance matrix
Omega=[1,0,0;0,1,0;0,0,1];

%Options
PERIOD=4;
HORIZON=60;

linlagex;
