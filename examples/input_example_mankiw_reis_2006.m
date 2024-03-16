%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%input_example_mankiw_reis_2006.m
%
%This is the model
%from: 
%Mankiw, N. G., and R. Reis (2006):"Pervasive Stickiness",
%American Economic Review, Papers and Proceedings, 96 (2), 164-169.
%
%
%THIS VERSION: 1.1 December 9, 2009: 
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('..\linlagex')

%Paramters
rho=0.92;
sigma_a=0.0085;
sigma_eps=0.0036;

lambda=1/4;
delta=1/4;
omega=1/4;


phi=4;
theta=1;
beta=2/3;
v=20;
gamma=10;

phi_pi=1.24;
phi_y=0.33;

%Steady States

%Paramter Abbreviations
eta_4a=((1+phi)*theta)/(phi*theta*(1-beta)+theta+beta*phi);

eta_5p=-(beta+v*(1-beta)*(1-lambda))/(beta+v*(1-beta));
eta_5w=lambda*beta/(beta+v*(1-beta));
eta_5y=lambda*(1-beta)/(beta+v*(1-beta));
eta_5a=-lambda/(beta+v*(1-beta));

eta_5p_it=v*(1-beta)*lambda*(1-lambda)/(beta+v*(1-beta));
eta_5w_it=beta*lambda*(1-lambda)/(beta+v*(1-beta));
eta_5y_it=(1-beta)*lambda*(1-lambda)/(beta+v*(1-beta));
eta_5a_it=-lambda*(1-lambda)/(beta+v*(1-beta));

eta_6yn=delta;
eta_6R=-theta*delta;

eta_6yn_it=delta*(1-delta);
eta_6R_it=-delta*(1-delta)*theta;

eta_7p=omega*phi/(gamma+phi);
eta_7w=-(gamma*(1-omega)+phi)/(gamma+phi);
eta_7y=omega*(1/beta)/(gamma+phi);
eta_7yn=omega*(phi/theta)/(gamma+phi);
eta_7a=-omega*(1/beta)/(gamma+phi);
eta_7R=-omega*phi/(gamma+phi);

eta_7p_it=eta_7p*(1-omega);
eta_7w_it=omega*(1-omega)*gamma/(gamma+phi);
eta_7y_it=eta_7y*(1-omega);
eta_7yn_it=eta_7yn*(1-omega);
eta_7a_it=eta_7a*(1-omega);
eta_7R_it=eta_7R*(1-omega);
%Variable names
ENDOGENOUS_VARIABLE_NAMES={'Long-Run Real Rate'
'Real Rate'
'Nominal Rate'
'Inflation'
'Output'
'Natural Output'
'Price Level'
'Nominal Wage'};
EXOGENOUS_VARIABLE_NAMES={'Technology'
'Interest Rate Shock'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0=sum A_iE_t-i(x_t+1)+sum B_iE_t-i(x_t)+sum C_iE_t-i(x_t-1)+sum
% F_iE_t-i(z_t+1)+sumG_iE_t-i(z_t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matrices
%       R   r   i   pi  y   y_r p   w 
A_0=[   1   0   0   0   0   0   0   0 
        0   0   0   1   0   0   0   0 
        0   0   0   0   0   0   0   0  
        0   0   0   0   0   0   0   0  
        0   0   0   0   0   0   0   0 
        0   0   0   0   0   0   0   0  
        0   0   0   0   0   0   0   0  
        0   0   0   0   0   0   0   0  ];
%       R   r   i   pi  y   y_r p   w  
B_0=[   -1  1   0   0   0   0   0   0
        0   1   -1  0   0   0   0   0  
        0   0   -1  phi_pi  phi_y -phi_y 0  0  
        0   0   0   0   0   -1  0   0   
        0   0   0   0   eta_5y  0   eta_5p  eta_5w  
        eta_6R   0   0   0   -1  eta_6yn 0   0   
        eta_7R  0   0   0   eta_7y  eta_7yn eta_7p  eta_7w  
        0   0   0   -1  0   0   1   0 ];
%       R   r   i   pi  y   y_r p   w 
C_0=[   0   0   0   0   0   0   0   0   
        0   0   0   0   0   0   0   0   
        0   0   0   0   0   0   0   0   
        0   0   0   0   0   0   0   0   
        0   0   0   0   0   0   0   0   
        0   0   0   0   0   0   0   0   
        0   0   0   0   0   0   0   0   
        0   0   0   0   0   0   -1   0  ];
%       a   eps
F_0=[	zeros(8,2)];
%       a   eps
G_0=[   0   0
        0   0
        0   1
        eta_4a  0
        eta_5a  0
        0   0
        eta_7a  0
        0   0];


%name of the counter used in subsequent matrices
it_name='J_1';
%maximum value of counter: either 'infinity' or a scalar (without
%quotation marks)
it_max_value='infinity';


%counter matrices

A_j='[zeros(8,8)]';
%        R   r   i   pi  y   y_r p   w
B_j=['[zeros(4,8);'...
    '0 0   0   0 eta_5y_it*(1-lambda)^(J_1-1) 0 eta_5p_it*(1-lambda)^(J_1-1) eta_5w_it*(1-lambda)^(J_1-1);'...
    'eta_6R_it*(1-delta)^(J_1-1) 0 0 0 0 eta_6yn_it*(1-delta)^(J_1-1) 0 0;'...
    'eta_7R_it*(1-omega)^(J_1-1) 0 0 0 eta_7y_it*(1-omega)^(J_1-1) eta_7yn_it*(1-omega)^(J_1-1) eta_7p_it*(1-omega)^(J_1-1) eta_7w_it*(1-omega)^(J_1-1);'...
    'zeros(1,8)]'];
C_j='[zeros(8,8)]';
F_j='[zeros(4,2);  eta_5a_it*(1-lambda)^(J_1-1) 0; 0 0; eta_7a_it*(1-omega)^(J_1-1) 0;0 0 ]';
G_j='[zeros(8,2)]';

%limit matrices
%If a counter goes to infinity, please enter the limiting summed matrices
clear A_inf B_inf C_inf F_inf G_inf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define exogenous process as a VAR(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AR-coefficients
N=[1,0;0,rho];
%Covariance matrix
Omega=[sigma_a^2,0;0,sigma_eps^2];

%Options
PERIOD=4;


linlagex;
