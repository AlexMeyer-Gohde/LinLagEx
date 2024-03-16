function [mean_correlations mean_impulses credible_corrs credible_imps median_corrs median_imps data_correlations]=bayesian_impulses_and_correlations(mixed_results, data, model_name, selection_matrix, ENDOGENOUS_VARIABLE_NAMES, EXOGENOUS_VARIABLE_NAMES, CORRELATION_HORIZON, IMPULSE_HORIZON, PERIOD, use_how_many,credible_set_percentile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%bayesian_impulses_and_correlations.m
%
%This program:
%               calculates and plots crosscorrelations from the posterior
%               distribution and the data and plots impulse responses from
%               the posterior distribution
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
%               CORRELATION_HORIZON: horizon of cross correlations 
%                   (i.e. corr(x_{t+j},y_t), where j is horizon)
%               IMPULSE_HORIZON: horizon out ot which impulse responses are
%                   to be calculated 
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
%               median_corrs: cell-by-cell median of crosscorrelations from posterior draws 
%               median_imps: ditto for impulses
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
[obs var]=size(data);
use_results_index=1:floor(length(mixed_results)/use_how_many):length(mixed_results);
use_results_index=use_results_index(1:use_how_many);
use_results=mixed_results(use_results_index,:);

total_correlations=cell(length(selection_matrix),length(selection_matrix));
total_impulses=cell(length(selection_matrix),length(selection_matrix));
for hh=1:use_how_many
    parameters=use_results(hh,:);
    [correlations impulses]=get_correlations(model_name,parameters,selection_matrix,CORRELATION_HORIZON,IMPULSE_HORIZON);
    for jjj=1:length(selection_matrix)
        for iii=1:length(selection_matrix)
            total_correlations{iii,jjj}=[total_correlations{iii,jjj}; correlations{iii,jjj}];
            total_impulses{iii,jjj}=[total_impulses{iii,jjj}; impulses{iii,jjj}];
        end
    end
end
credible_corrs=cell(length(selection_matrix),length(selection_matrix));
credible_imps=cell(length(selection_matrix),length(selection_matrix));
median_corrs=cell(length(selection_matrix),length(selection_matrix));
median_imps=cell(length(selection_matrix),length(selection_matrix));
for jjj=1:length(selection_matrix)
        for iii=1:length(selection_matrix)
            for j=1:2*CORRELATION_HORIZON+1
                temp= sort(total_correlations{iii,jjj}(:,j));
                credible_corrs{iii,jjj}(1,j)=temp(round(credible_set_percentile*length(total_correlations{iii,jjj})+1));
                credible_corrs{iii,jjj}(2,j)=temp(round((1-credible_set_percentile)*length(total_correlations{iii,jjj})));
                median_corrs{iii,jjj}(1,j)=median(total_correlations{iii,jjj}(:,j));
            end
            for j=1:IMPULSE_HORIZON
                temp= sort(total_impulses{iii,jjj}(:,j));
                credible_imps{iii,jjj}(1,j)=temp(round(credible_set_percentile*length(total_impulses{iii,jjj})+1));
                credible_imps{iii,jjj}(2,j)=temp(round((1-credible_set_percentile)*length(total_impulses{iii,jjj})));
                median_imps{iii,jjj}(1,j)=median(total_impulses{iii,jjj}(:,j));
            end
        end
end
[mean_correlations mean_impulses]=get_correlations(model_name,mean(mixed_results),selection_matrix,CORRELATION_HORIZON,IMPULSE_HORIZON);




for j=0:CORRELATION_HORIZON
Gamma_1_data{j+1}=data(1:obs-j,:)'*data(1+j:obs,:);
end
for j=1:2*length(Gamma_1_data)-1
for i=1:var
if j-length(Gamma_1_data)<0
Xcorr_data{i}(:,j)=Gamma_1_data{abs(j-(length(Gamma_1_data)+1))}(:,i)./((Gamma_1_data{1}(i,i).^(1/2))*(diag(Gamma_1_data{1}).^(1/2)));
else
Xcorr_data{i}(:,j)=Gamma_1_data{abs(j-(length(Gamma_1_data))+1)}(i,:)'./((Gamma_1_data{1}(i,i).^(1/2))*(diag(Gamma_1_data{1}).^(1/2)));
end
end
end
for jjj=1:length(selection_matrix)
for iii=1:length(selection_matrix)
    data_correlations{iii,jjj}=Xcorr_data{selection_matrix(jjj)}(selection_matrix(iii),1:2*CORRELATION_HORIZON+1);
end
end

for iii=1:length(ENDOGENOUS_VARIABLE_NAMES)
for jjj=1:length(ENDOGENOUS_VARIABLE_NAMES)
figure;
plot(-CORRELATION_HORIZON:CORRELATION_HORIZON, data_correlations{iii,jjj},'-k',-CORRELATION_HORIZON:CORRELATION_HORIZON, mean_correlations{iii,jjj},':k','LineWidth',4)
eval(sprintf('title(''Cross-Correlations of %s at t+j with %s at t'')',ENDOGENOUS_VARIABLE_NAMES{iii},ENDOGENOUS_VARIABLE_NAMES{jjj}))
legend({'Data','Posterior'},'Location','Best');
ylabel('Correlation Coefficient');
xlabel('j');
axis([-Inf,Inf,-1,1])
hold on
plot(-CORRELATION_HORIZON:CORRELATION_HORIZON, credible_corrs{iii,jjj}(1,:),':k',-CORRELATION_HORIZON:CORRELATION_HORIZON, credible_corrs{iii,jjj}(2,:),':k','LineWidth',2)
end
end


TIME=(-1*PERIOD:IMPULSE_HORIZON-1)/PERIOD;
for iii=1:length(ENDOGENOUS_VARIABLE_NAMES)
for jjj=1:length(EXOGENOUS_VARIABLE_NAMES)
figure;
plot(TIME,[zeros(1,PERIOD)  mean_impulses{iii,jjj}],':k','LineWidth',4)
eval(sprintf('title(''Impulse Responses of %s to a Unit Shock in %s'')',ENDOGENOUS_VARIABLE_NAMES{iii},EXOGENOUS_VARIABLE_NAMES{jjj}))
ylabel('% Deviations from Steady-State');
xlabel('Years since Shock Realization');
hold on
plot(TIME, 0*TIME,'k')
plot(TIME,[zeros(1,PERIOD) credible_imps{iii,jjj}(1,:)],':k',TIME,[zeros(1,PERIOD) credible_imps{iii,jjj}(2,:)],':k','LineWidth',2)
end
end