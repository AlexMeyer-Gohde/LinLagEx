%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%input_example_mankiw_reis_2002.m
%
%This is the model
%from: 
%Mankiw, N. G., and R. Reis (2002):"Sticky Information versus Sticky Prices: A Proposal to Replace the New Keynesian Phillips Curve",
%The Quarterly Journal of Economics, 117(4), 1295-1328.
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

addpath('..\linlagex')

%Paramters
alpha=0.1;
lambda=0.25;
rho_m=0.5;
sigma_m=0.007;



%Steady States



%Variable names
ENDOGENOUS_VARIABLE_NAMES={'Inflation'
'Output'
'Output Growth'};
EXOGENOUS_VARIABLE_NAMES={'Money Growth'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0=sum A_iE_t-i(x_t+1)+sum B_iE_t-i(x_t)+sum C_iE_t-i(x_t-1)+sum
% F_iE_t-i(z_t+1)+sumG_iE_t-i(z_t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matrices
%       pi  y   dy
A_0=[   zeros(3,3)];
%       pi  y   dy
B_0=[   1   0   1
        -1 alpha*lambda/(1-lambda) 0
        0   -1  1];
%       pi  y   dy
C_0=[   0   0   0
        0   0   0
        0   1   0];
%       dm
F_0=[	zeros(3,1)];
%       dm
G_0=[   -1*(-sigma_m)
        zeros(2,1)];


%name of the counter used in subsequent matrices
it_name='J_1';
%maximum value of counter: either 'infinity' or a scalar (without
%quotation marks)
it_max_value='infinity';


%counter matrices

A_j='[zeros(3,3)]';
%                   pi              y  dy
B_j='[zeros(1,3);lambda*(1-lambda)^(J_1-1) 0  alpha*lambda*(1-lambda)^(J_1-1); zeros(1,3)]';
C_j='[zeros(3,3)]';
F_j='[zeros(3,1)]';
G_j='[zeros(3,1)]';

%limit matrices
%If a counter goes to infinity, please enter the limiting summed matrices
A_inf=[A_0];
B_inf=[   1   0   1
        0 alpha*lambda/(1-lambda) alpha
        0   -1  1];
C_inf=[C_0];
F_inf=[F_0];
G_inf=[G_0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define exogenous process as a VAR(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AR-coefficients
%       da  dm
N=[rho_m];
%Covariance matrix
Omega=[sigma_m^2];

%Options
PERIOD=4;

linlagex;
