function [mean_decomps median_decomps credible_decomps]=bayesian_variance_decomposition(mixed_results, model_name, selection_matrix, ENDOGENOUS_VARIABLE_NAMES, EXOGENOUS_VARIABLE_NAMES, decomp_horizon, PERIOD, use_how_many,credible_set_percentile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%bayesian_variance_decompositions.m
%
%This program:
%               calculates variance decompositions at user defined forecast 
%               horizons (plus the infinite horizon) from the posterior distribution
%
%Inputs:        mixed_results: a matrix containing the draws from the
%                   posterior (each row corresponds to a draws and each column
%                   corresponds to a parameter)
%               data: the data (row correspond to time in ascending order and 
%                   columns corrspond to observables) 
%               model_name: a string containing the name of the model file 
%               selection_matrix:: a vector containing an integer at the
%                   i'th entry that denotes which endogenous variable in
%                   the model file corresponds to the i'th observable
%               ENDOGENOUS_VARIABLE_NAMES: a cell vector of strings, 
%                   row i contains the name of the observable in column i
%                   of the data matrix
%               EXOGENOUS_VARIABLE_NAMES: a cell vector of strings
%                   containing the name of the exogenous shocks
%               decomp_horizon: forecast horizons for variance decomps
%               PERIOD: how many periods per year 
%               use_how_many: how many of the draws to use. I.e., if there are 100 draws 
%                   and use_how_many=10, then draw 1, 11, 21,.., 91 will be used
%               credible_set_percentile: an number between 0 and 1
%                   corresponding to the percentiles to be used for the
%                   credible set. I.e. if credible_set_percentile=0.05,
%                   then the 5% and 95% bounds will be calculated
%              
%
%Output:        mean_correlations: crosscorrelations at the mean of mixed_results
%               mean_impulses: impulses at ditto 
%               credible_corrs: crosscorrelations at the percentiles
%                   determined bycredible_set_percentile
%               credible_imps: impulses at ditto 
%               data_correlations: crosscorrelations of the data
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
model_name=str2func(model_name);


IMPULSE_HORIZON=max(decomp_horizon);



use_results_index=1:floor(length(mixed_results)/use_how_many):length(mixed_results);
use_results_index=use_results_index(1:use_how_many);
use_results=mixed_results(use_results_index,:);

total_decomps=cell(length(selection_matrix),length(decomp_horizon)+1);


for hh=1:use_how_many
    parameters=use_results(hh,:);
    [decomps]=get_decomps(model_name,parameters,selection_matrix,decomp_horizon,IMPULSE_HORIZON);
    for iii=1:length(selection_matrix)
        for jjj=1:length(decomp_horizon)+1
            total_decomps{iii,jjj}=[total_decomps{iii,jjj}; decomps{iii,jjj}];
        end
    end
end



credible_decomps=cell(length(selection_matrix),length(decomp_horizon)+1);
median_decomps=cell(length(selection_matrix),length(decomp_horizon)+1);
for jjj=1:length(decomp_horizon)+1
        for iii=1:length(selection_matrix)
            for j=1:length(EXOGENOUS_VARIABLE_NAMES)
                temp= sort(total_decomps{iii,jjj}(:,j));
                credible_decomps{iii,jjj}(1,j)=temp(round(credible_set_percentile*length(total_decomps{iii,jjj})+1));
                credible_decomps{iii,jjj}(2,j)=temp(round((1-credible_set_percentile)*length(total_decomps{iii,jjj})));
                median_decomps{iii,jjj}(1,j)=median(total_decomps{iii,jjj}(:,j));
            end
            end
end
[mean_decomps]=get_decomps(model_name,mean(mixed_results),selection_matrix,decomp_horizon,IMPULSE_HORIZON);



