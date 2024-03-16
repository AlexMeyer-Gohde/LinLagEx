function [decomps]=get_decomps(model_name,parameters,selection_matrix,decomp_horizon,IMPULSE_HORIZON);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%get_decomps.m
%
%gets variance decomps for specified (plus infinite) forecast horizons.
%
%THIS VERSION: 1.1 December 9, 2009
%
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DISPLAY_STD=0;
PLOT_CORRELATION=0;
IMPULSE_SELECT=selection_matrix;
HORIZON=IMPULSE_HORIZON;
model_name();
 for iii=1:length(selection_matrix)
        for jjj=1:length(decomp_horizon)
All_Var=diag((chol(Omega)*Resp_mat(([0:length(EXOGENOUS_VARIABLE_NAMES)-1])*(num_endog+num_exog)+IMPULSE_SELECT(iii),1:decomp_horizon(jjj)))*(chol(Omega)*Resp_mat(([0:length(EXOGENOUS_VARIABLE_NAMES)-1])*(num_endog+num_exog)+IMPULSE_SELECT(iii),1:decomp_horizon(jjj)))');
Total_Var=sum(All_Var);
decomps{iii,jjj}=All_Var'/Total_Var;
end
end
Omega_store=Omega;
spectral
Total_Var=diag(Gamma_1_spec{length(s_Y)/2+1});
Total_Var=Total_Var(selection_matrix);
Decomp_Var=[];
for jjj=1:num_exog
Omega=zeros(num_exog,num_exog);
Omega(jjj,jjj)=Omega_store(jjj,jjj);
spectral
All_Var=diag(Gamma_1_spec{length(s_Y)/2+1});
Decomp_Var=[Decomp_Var All_Var(selection_matrix)./Total_Var];
end
for iii=1:length(selection_matrix)
decomps{iii,length(decomp_horizon)+1}=Decomp_Var(iii,:);
end