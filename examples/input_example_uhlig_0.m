%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%input_example_uhlig_0.m
%
%This is exampl0.m from Prof. Uhlig's Toolkit
%
% VERSION 2.0, MARCH 1997, COPYRIGHT H. UHLIG.
% EXAMPL0.M:
% Solving the stochastic neoclassical growth model with the "toolkit"
% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.
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

% Setting parameters:

Z_bar     = 1;    % Normalization
rho       = .36;  % Capital share
delta     = .025; % Depreciation rate for capital
R_bar     = 1.01; % One percent real interest per quarter
eta       = 1.0;  % constant of relative risk aversion = 1/(coeff. of intertemporal substitution)
psi       = .95;  % autocorrelation of technology shock
sigma_eps = .712; % Standard deviation of technology shock.  Units: Percent.

% Calculating the steady state:

betta   = 1.0/R_bar;  % Discount factor beta
K_bar   = ((rho*Z_bar)/(R_bar - 1 + delta))^(1.0/(1 - rho));
Y_bar   = Z_bar*K_bar^rho;
C_bar   = Y_bar - delta*K_bar;

%Variable names
ENDOGENOUS_VARIABLE_NAMES={'capital'
            'consumption'
            'return'
            'output'};
EXOGENOUS_VARIABLE_NAMES={'technology'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0=sum A_iE_t-i(x_t+1)+sum B_iE_t-i(x_t)+sum C_iE_t-i(x_t-1)+sum
% F_iE_t-i(z_t+1)+sumG_iE_t-i(z_t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matrices

A_0=[0 0 0 0
    0 0 0 0
    0 0 0 0
    0 -eta 1 0];
B_0=[0 0 -1 0
     -K_bar/C_bar -1 0 0
      0 0 0 -1
      0 eta 0 0];
C_0=[-(1-betta*(1-delta))*(1-rho) 0 0 0
    K_bar/(betta*C_bar) 0 0 0
    rho 0 0 0
    0 0 0 0];
F_0=[zeros(4,1)];
G_0=[(1-betta*(1-delta))
      (1+delta*K_bar/C_bar)
    1
    0];


%name of the counter used in subsequent matrices
%it_name='J_1';
%maximum value of counter: either 'infinity' or a scalar (without
%quotation marks)
it_max_value=0;


%counter matrices



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define exogenous process as a VAR(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AR-coefficients
N=[psi];
%Covariance matrix
Omega=[sigma_eps^2];

%Options
PERIOD=4;


linlagex;
