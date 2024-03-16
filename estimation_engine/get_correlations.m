function [correlations impulses]=get_correlations(model_name,parameters,selection_matrix,CORRELATION_HORIZON,IMPULSE_HORIZON);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%get_correlations.m
%
%This program exports crosscorrelations and impulses for use in
%bayesian_impulses_and_correlations
%
%Inputs:        model_name: a function that executes the model file
%               parameters: a vector of parameters to be passed to the
%                   model
%               selection_matrix: a vector containing an integer at the
%                   i'th entry that denotes which endogenous variable in
%                   the model file corresponds to the i'th observable
%               CORRELATION_HORIZON: horizon of cross correlations (i.e. corr(x_{t+j},y_t), where j is horizon)
%               IMPULSE_HORIZON: horizon out ot which impulse responses are
%                   to be calculated
%
%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CORRELATION_SELECT=selection_matrix;
IMPULSE_SELECT=selection_matrix;
HORIZON=IMPULSE_HORIZON;
model_name();
for jjj=1:length(selection_matrix)
for iii=1:length(selection_matrix)
    correlations{iii,jjj}=Xcorr_spec{CORRELATION_SELECT(jjj)}(CORRELATION_SELECT(iii),grid_size/2+1-CORRELATION_HORIZON:grid_size/2+1+CORRELATION_HORIZON);
    impulses{iii,jjj}=Resp_mat((jjj-1)*(num_endog+num_exog)+IMPULSE_SELECT(iii),:);

end
end
