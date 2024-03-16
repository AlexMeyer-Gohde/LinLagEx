%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%posterior_mode_sp.m
%
%This program:
%               is a script to get the posterior mode approximated by
%               numerically maximizing the log posterior likelihood
%               function for the sticky-price model
%
%       Simply change to the directory where this file is located, type
%       posterior_mode_sp and the results will be saved in the file
%       posterior_mode_sp.mat
%
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

load sw_2005_full;  
x=data(:,2);
hp_filter; %HP Filter the GDP data
data(:,2)=xhp;
data(:,1)=data(:,1)-mean(data(:,1));    %De-mean inflation
data(:,3)=data(:,3)-mean(data(:,3));    %De-mean interest rates
[T K]=size(data);

model_name='input_example_me_sp';

Pmean=[0.5 0.25      1   1.5   0.125 0.5  0.5 0.5        0.1 0.1 0.1]; 
Pstd= [0.2 0.05  0.75   0.125  0.05  0.2   0.2  0.2 2      2 2];
Plower=[0   0      0      1     0       0   0   0  0   0         0 ];
Pupper=[1   1    1/0     1/0    1/0    1   1   1   1/0  1/0    1/0];
Ptype={'beta', 'beta', 'norm', 'norm', 'norm', 'beta', 'beta', 'beta', 'invgamma2', 'invgamma2', 'invgamma2'};
Pparameter=[Pmean;Pstd;Plower;Pupper];



Q_1=zeros(3,3);
Q_1(1,1)=1;
Q_1(2,2)=1;
Q_1(3,3)=1;
Q{1}=Q_1;
Q{2}=[];

save_name='posterior_mode_sp.mat';

x0=Pmean;



Data=data';
Data=Data(:);
addpath('..\estimation_engine')
addpath('..\estimation_engine\sims_opt')
addpath('..\linlagex')
mode_estimates(Pparameter,Ptype,Q,Data,T,K,model_name,x0,save_name)