%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%input_example_trabandt_2007.m
%
%This is the sticky information model
%from: 
%Trabandt, M. (2007): “Sticky Information vs. Sticky Prices: A Horse Race
%in a DSGE Framework,” Sveriges Riksbank Working Paper Sieres 209.
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

%Parameters
lambda=0.75;

rho_z=0.95;
sigma_eps_z=0.0071;
rho_g=0.95;
sigma_eps_g=0.006;
rho_zi=0.50;
sigma_eps_zi=0.0080;

beta=0.99;
sigma=2;
phi=1.5;
v=2;
alpha=2/3;
theta=6;

%Steady States
g_bar=0.3;
zT_bar=1;
R_bar=1/beta;
y_bar=zT_bar;
c_bar=y_bar-g_bar;
sc=c_bar/y_bar;
w=phi/alpha+1/alpha-1;
xi=(w+sigma/sc)/(1+theta*w);
a_1=sc/sigma;
muz=sigma*(1+w)*(rho_z-1)/(sc*w+sigma);
a_3=sc*muz/sigma;
mug=(sigma*(rho_g-1)/sc)*((sigma*(1-sc)/(sc*w+sigma))+sc-1);
a_2=sc*mug/sigma;
psi=(w+sigma/sc)/(1+w);
gamma_g=(sigma*(1-sc)/(sc*v))*(1-(sigma/sc)/(w+sigma/sc));


%Variable names
ENDOGENOUS_VARIABLE_NAMES={'Real Money Supply'
'Inflation'
'Natural Real Interest Rate'
'Output Gap'
'Nominal Interest Rate'
'Change in Output Gap'};
EXOGENOUS_VARIABLE_NAMES={'Technology'
'Government Spending'
'Money Supply'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0=sum A_iE_t-i(x_t+1)+sum B_iE_t-i(x_t)+sum C_iE_t-i(x_t-1)+sum
% F_iE_t-i(z_t+1)+sumG_iE_t-i(z_t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matrices
%       m   pi  rr  x   R   Dx
A_0=[   0   0   0   0   0   0
        0   0   0   0   0   0
        0   0   0   0   0   0
        0   sc/sigma    0   1   0   0
        0   0   0   0   0   0
        0   0   0   0   0   0];
%       m   pi  rr  x   R   Dx
B_0=[   -1  -1  0   0   0   0
        0   0   -1  0   0   0
        -1  0   0   sigma/(sc*v)    -1/(v*(R_bar-1))    0
        0   0   sc/sigma    -1  -sc/sigma   0
        0   -1  0   xi*(1-lambda)/lambda    0   0
        0   0   0   1   0   -1];
%       m   pi  rr  x   R   Dx
C_0=[   1   0   0   0   0   0
        0   0   0   0   0   0
        0   0   0   0   0   0
        0   0   0   0   0   0
        0   0   0   0   0   0
        0   0   0   -1  0   0];
%       z   g   zi
F_0=[	zeros(6,3)];
%       z   g   zi
G_0=[   0   0   1
        muz mug 0
        sigma/(sc*v*psi)    -gamma_g    0
        0   0   0
        0   0   0
        0   0   0];


%name of the counter used in subsequent matrices
it_name='J_1';
%maximum value of counter: either 'infinity' or a scalar (without
%quotation marks)
it_max_value='infinity';


%counter matrices

A_j='[zeros(6,6)]';
%       m   pi  rr  x   R   Dx
B_j='[zeros(4,6); 0 (1-lambda)*lambda^(J_1-1)  0 0 0 xi*(1-lambda)*lambda^(J_1-1); zeros(1,6)]';
C_j='[zeros(6,6)]';
F_j='[zeros(6,3)]';
G_j='[zeros(6,3)]';

%limit matrices
%If a counter goes to infinity, please enter the limiting summed matrices
A_inf=[A_0];
B_inf=[-1  -1  0   0   0   0
        0   0   -1  0   0   0
        -1  0   0   sigma/(sc*v)    -1/(v*(R_bar-1))    0
        0   0   sc/sigma    -1  -sc/sigma   0
        0   0  0   xi*(1-lambda)/lambda    0   xi
        0   0   0   1   0   -1];
C_inf=[C_0];
F_inf=[F_0];
G_inf=[G_0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define exogenous process as a VAR(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AR-coefficients
N=[rho_z,0,0;0,rho_g,0;0,0,rho_zi];
%Covariance matrix
Omega=[sigma_eps_z^2,0,0; 0,sigma_eps_g^2,0;0,0,sigma_eps_zi^2];

%Options
PERIOD=4;


linlagex;
