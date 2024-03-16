%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%anticipated_impulse.m
%
%
%As the solution was formulated as an infinite vector MA,
%impulse responses to anticipated/pre-announced innovations can be solved for
% a la Taylor (1986)
%
%As in the program impulse_response.m, the basic programming structure of:
%The following text and subsequent basic structure of the program has been taken/adapted from: 
% VERSION 4.0, November 2002, COPYRIGHT H. UHLIG.
% IMPRESP.M
%has been used here.
%
%THIS VERSION: 1.1 December 9, 2009
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('ANTICIPATED_HORIZON','var')==0
ANTICIPATED_HORIZON=40;
end
if exist('ANTICIPATED_PERIOD','var')==0
ANTICIPATED_PERIOD=8;
end


if ANTICIPATED_PERIOD<=max_it
DIM_DUM=max_it;
RHS=zeros((max_it+1)*num_eqs,num_exog);
F=F_0;
G=G_0;
for j=1:max_it-1
F=F+F_j(j);
G=G+G_j(j);
if j>=ANTICIPATED_PERIOD
RHS(1+j*num_eqs:num_eqs+j*num_eqs,1:num_exog)=-(F*(N^(j+1-ANTICIPATED_PERIOD))+G*(N^(j-ANTICIPATED_PERIOD)));
end
end
RHS(1+max_it*num_eqs:(max_it+1)*num_eqs,1:num_exog)=BETA_ZS*N^(max_it-ANTICIPATED_PERIOD);

else
DIM_DUM=ANTICIPATED_PERIOD;
if max_it==0
LHS=[sparse((ANTICIPATED_PERIOD+1)*num_eqs,(ANTICIPATED_PERIOD+1)*num_endog)];
C=C_0;
B=B_0;
A=A_0;
else
LHS=[LHS(1:max_it*num_eqs,:),sparse((max_it)*num_eqs,(ANTICIPATED_PERIOD-max_it)*num_endog);...
sparse((ANTICIPATED_PERIOD-max_it+1)*num_eqs,(ANTICIPATED_PERIOD+1)*num_endog)];
end
RHS=zeros((ANTICIPATED_PERIOD+1)*num_eqs,num_exog);
for j=max_it:ANTICIPATED_PERIOD-1
if j>0
LHS(1+j*num_eqs:num_eqs+j*num_eqs,1+(j-1)*num_endog:(j+2)*num_endog)=[C,B,A];
else
LHS(1:num_eqs,1:2*num_endog)=[B_0,A_0];
end
end
LHS(1+ANTICIPATED_PERIOD*num_eqs:(ANTICIPATED_PERIOD+1)*num_eqs,1+(ANTICIPATED_PERIOD-1)*num_endog:(ANTICIPATED_PERIOD+1)*num_endog)=[-ALPHA_ZS,eye(num_endog,num_endog)];
RHS(1+ANTICIPATED_PERIOD*num_eqs:(ANTICIPATED_PERIOD+1)*num_eqs,1:num_exog)=BETA_ZS;
end
ANTICIPATED_IMPULSE_VECTOR=LHS\RHS;
ANTICIPATED_IMPULSE_VECTOR=full(ANTICIPATED_IMPULSE_VECTOR);



ANTICIPATED_RESP_MAT = [];
ANTICIPATED_RESPONSE= [];
EXOG_RESPONSE_TEMP=[];
EXOG_RESPONSE=[];
for j=1:max(ANTICIPATED_HORIZON,DIM_DUM+1)
if j>ANTICIPATED_PERIOD
EXOG_RESPONSE_TEMP=N^(j-1-ANTICIPATED_PERIOD);
else
EXOG_RESPONSE_TEMP=zeros(size(N));
end
EXOG_RESPONSE=[EXOG_RESPONSE;EXOG_RESPONSE_TEMP];
end



if ANTICIPATED_HORIZON<=DIM_DUM+1
    for j=1:num_exog
        ANTICIPATED_RESPONSE=reshape(ANTICIPATED_IMPULSE_VECTOR(:,j),[num_endog,DIM_DUM+1]);
        EXOG_RESPONSE_TEMP=reshape(EXOG_RESPONSE(1:num_exog*(DIM_DUM+1),j),[num_exog,DIM_DUM+1]);
        ANTICIPATED_RESP_MAT=[ANTICIPATED_RESP_MAT;ANTICIPATED_RESPONSE;EXOG_RESPONSE_TEMP];
    end
    ANTICIPATED_RESP_MAT=ANTICIPATED_RESP_MAT(:,1:ANTICIPATED_HORIZON);
elseif ANTICIPATED_HORIZON>DIM_DUM+1
        ANTICIPATED_RESPONSE=[ANTICIPATED_IMPULSE_VECTOR];
    
        ANTICIPATED_IMPULSE_VECTOR_LAST=ANTICIPATED_IMPULSE_VECTOR(DIM_DUM*num_endog+1:(DIM_DUM+1)*num_endog,:);
    
    for j=DIM_DUM+2:ANTICIPATED_HORIZON
        ANTICIPATED_IMPULSE_VECTOR_LAST=ALPHA_ZS*ANTICIPATED_IMPULSE_VECTOR_LAST+BETA_ZS*N^(j-(DIM_DUM+1)-min(max_it,ANTICIPATED_PERIOD)+max_it);
        ANTICIPATED_RESPONSE=[ANTICIPATED_RESPONSE;ANTICIPATED_IMPULSE_VECTOR_LAST];
    end
    for j=1:num_exog
        ANTICIPATED_RESPONSE_temp=reshape(ANTICIPATED_RESPONSE(:,j),[num_endog,ANTICIPATED_HORIZON]);
        EXOG_RESPONSE_TEMP=reshape(EXOG_RESPONSE(:,j),[num_exog,ANTICIPATED_HORIZON]);
        ANTICIPATED_RESP_MAT=[ANTICIPATED_RESP_MAT;ANTICIPATED_RESPONSE_temp;EXOG_RESPONSE_TEMP];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plots for anticipated impulse responses %%%%%%%%%
if exist('PLOT_ANTICIPATED_IMPULSE','var')==0 || PLOT_IMPULSE==1 %If either the option has not been set or if the option PLOT_ANTICIPATED_IMPULSE has been set to 1
% If the user didn't specify a subset of variables to be displayed, all
% variables are displayed
if exist('ANTICIPATED_SELECT','var')==0
ANTICIPATED_SELECT=1:(num_endog+num_exog);
end

% If the user didn't specify the number of periods in a year, quarters are
% assumed
if exist('PERIOD','var')==0
PERIOD=4;
end


TIME=(-1*PERIOD:ANTICIPATED_HORIZON-1)/PERIOD;
for j=1:num_exog
    clf('reset')
    set(gcf,'DefaultAxesColorOrder',[hsv(length(ANTICIPATED_SELECT)-1);0.25, 0.25,0.25])
    IMPULSE=plot(TIME,[zeros(length(ANTICIPATED_SELECT),PERIOD) ANTICIPATED_RESP_MAT((j-1)*(num_endog+num_exog)+ANTICIPATED_SELECT,:)]',':.','MarkerEdgeColor','auto','MarkerSize',8);
    eval(sprintf('title({''Impulse Responses to a Shock in %s''; ''Anticipated/Pre-Announced %d Period(s) in Advance''})',EXOGENOUS_VARIABLE_NAMES{j},ANTICIPATED_PERIOD))
    if (ANTICIPATED_PERIOD)/PERIOD<=(ANTICIPATED_HORIZON-1)/PERIOD
        hold on
        plot((ANTICIPATED_PERIOD)/PERIOD+zeros(1,11),min(ylim):(max(ylim)-min(ylim))/10:max(ylim),'k:');
        hold off
        legend([VARIABLE_NAMES(ANTICIPATED_SELECT,1);'Date of Shock Realization'],'Location','Best');
    else
        legend([VARIABLE_NAMES(ANTICIPATED_SELECT,1)],'Location','Best');
    end
    ylabel('% Deviations from Steady-State');
    xlabel('Years since Announcement');
    hold on
    plot(TIME, 0*TIME,'k')
    hold off
    pause;
    clf('reset')
    for i=1:length(VARIABLE_NAMES)
        subplot(ceil(length(VARIABLE_NAMES)^(1/2)),round(length(VARIABLE_NAMES)^(1/2)),i);
        plot(TIME, [zeros(1,PERIOD) ANTICIPATED_RESP_MAT((j-1)*(num_endog+num_exog)+i,:)],'k:.','MarkerEdgeColor','auto','MarkerSize',8);
        
        if (ANTICIPATED_PERIOD)/PERIOD<=(ANTICIPATED_HORIZON-1)/PERIOD
            hold on
            plot((ANTICIPATED_PERIOD)/PERIOD+zeros(1,11),min(ylim):(max(ylim)-min(ylim))/10:max(ylim),'k:');
            hold off
            legend([VARIABLE_NAMES(i,1);'Date of Shock Realization'],'Location','Best');
        else
            legend(VARIABLE_NAMES(i,1),'Location','Best');
        end
        hold on
        plot(TIME, 0*TIME,'k')
        hold off
        ylabel('% Deviations from Steady-State');
        xlabel('Years since Announcement');
        eval(sprintf('title({''Impulse Responses to a Shock in %s''; ''Anticipated/Pre-Announced %d Period(s) in Advance''})',EXOGENOUS_VARIABLE_NAMES{j},ANTICIPATED_PERIOD))
    end
    pause;
end
end