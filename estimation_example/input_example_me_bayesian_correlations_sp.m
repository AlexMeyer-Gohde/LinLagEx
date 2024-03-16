%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%input_example_me_bayesian_correlations_truncation_sp.m
%
%This program:
%               is the model file used for calculating posterior
%               correlations and impulses in the sticky-price model.
%
%           The only difference to
%           input_example_me_sp.m is in the options: we want to calculate
%           impulses and frequency domain moments, but without dsplaying
%           anything, as we'll be re-running this file for a multitude of
%           values in the parameter vector parameter.
%
%           Note the mapping from the vector "parameters" to the model
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


beta=0.99;


%Variable names
ENDOGENOUS_VARIABLE_NAMES={'Inflation'
'Output Gap'
'Nominal Interest Rate'};
EXOGENOUS_VARIABLE_NAMES={'IS'
'PC'
'MP'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0=sum A_iE_t-i(x_t+1)+sum B_iE_t-i(x_t)+sum C_iE_t-i(x_t-1)+sum
% F_iE_t-i(z_t+1)+sumG_iE_t-i(z_t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matrices

A_0=[beta 0 0  
    a_1 1 0  
    0 0 0  ];
B_0=[-1 (1-lambda)*(1-lambda*beta)*xi/lambda 0  
    0 -1 -a_1  
(1-phi_R)*phi_pi   (1-phi_R)*phi_y   -1  ];
C_0=[zeros(2,3) 
        0 0 phi_R  ];
F_0=[zeros(3,3)];
G_0=[0 (1-lambda)*(1-lambda*beta)/lambda 0
    1 0 0
0 0 1];


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
RUN_IMPULSE=1;
RUN_ANTICIPATED_IMPULSE=0;
RUN_SIMULATION=0;
RUN_SPECTRAL=1;
PLOT_IMPULSE=0;
PLOT_CORRELATION=0;
DISPLAY_STD=0;
solution_method='AIM';
linlagex;