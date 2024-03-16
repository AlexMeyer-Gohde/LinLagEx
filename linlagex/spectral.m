%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%spectral.m
%
%This program uses a recursive representation of the infinite MA solution
%form to calculate population moments.
%
%The code for the grid points, the HP transfer function and the basic
%programming structure for calculating moments via the population spectrum
%have been taken/adapted from:
% VERSION 4.1, MAY 2003, COPYRIGHT H. UHLIG.   
% MOMENTS.M calculates variances, covariances, and autocorrelation 
% tables, using frequency domain techniques with and without HP-Filtering.
%
%Moments are calculated to allow for
%a broader (w.r.t. H. UHLIG's toolkit) range of cross-correlations 
%to be calculated
%
%THIS VERSION: 1.1 December 9, 2009
%
%Fixed: sim references 9.1.2009, spectral calculations
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gamma_1_spec=[];
Gamma_1_spec_hp=[];
Xcorr_spec=[];
Xcorr_spec_hp=[];


if max(abs(eig(N)))>=1 || max(abs(eig(ALPHA_ZS)))>=1
disp('spectral.m')
disp('Sorry, the system is not covariance stationary')
disp('I cannot calculate population moments via')
disp('frequency domain methods')
disp('')
disp('You can still try to obtain HP-filtered sample moments')
disp('by running a simulation (RUN_SIMULATION=1)')
else
if exist('grid_size','var')==0
grid_size=64;
end
if exist('HP_LAMBDA','var')==0
if exist('PERIOD','var')==0
disp('You did not set PERIOD (number of periods per year)')
disp('I shall set HP_LAMBDA under the assumption that your model')
disp('is quarterly (four periods per year).')
PERIOD=4;
end
HP_LAMBDA=1600*(PERIOD/4)^4; %Ravn, Morten and Harald Uhlig (2002) "On Adjusting the HP-Filter for the Frequency of Observations"
end


grid_points=0:2*pi/grid_size:2*pi-2*pi/grid_size;


s_Y=cell(1,grid_size);
s_Y_hp=cell(1,grid_size);
s_yy=zeros((num_endog+num_exog)*(num_endog+num_exog),grid_size);
s_yy_hp=zeros((num_endog+num_exog)*(num_endog+num_exog),grid_size);
Gamma_1_spec=cell(1,length(s_Y)+1);
Gamma_1_spec_hp=cell(1,length(s_Y)+1);
Xcorr_spec=cell(1,num_endog+num_exog);
Xcorr_spec_hp=cell(1,num_endog+num_exog);

hp = 4*HP_LAMBDA*(1 - cos(grid_points)).^2 ./ (1 + 4*HP_LAMBDA*(1 - cos(grid_points)).^2);

%reshape block column vector [theta(0); theta(1); ...] into block row
%vector [theta(0), theta(1), ...]
THETA=reshape(MA_VECTOR(1:num_endog*max_it,:), [ num_endog max_it num_exog ] );
THETA=permute(THETA, [1 3 2]);
THETA=reshape(THETA,[num_endog num_exog*(max_it)]);


if max_it>0
THETA_I_MINUS_1=MA_VECTOR((max_it-1)*num_endog+1:(max_it)*num_endog,:);
else
THETA_I_MINUS_1=zeros(num_endog,num_exog);
end

for n = 1 : grid_size
e_plus=exp( ((-1)^(1/2))*grid_points(n));
E_PLUS=kron(e_plus.^(0:max_it-1),eye(num_exog,num_exog));
e_minus=exp( (-(-1)^(1/2))*grid_points(n));
E_MINUS=kron(e_minus.^(0:max_it-1),eye(num_exog,num_exog));
s_Y{n}=((2*pi)^(-1))*([THETA*conj(E_MINUS')+(inv(eye(num_endog)-ALPHA_ZS*e_minus)*...
        (BETA_ZS*(N^(max_it))*inv(eye(num_exog)-N*e_minus)+ALPHA_ZS*THETA_I_MINUS_1))*(e_minus^(max_it));...
        inv(eye(num_exog)-N*e_minus)]*Omega*[E_PLUS*THETA'+(e_plus^(max_it))*...
        (inv(eye(num_exog)-N'*e_plus)*(N^(max_it)')*BETA_ZS'+THETA_I_MINUS_1'*ALPHA_ZS')*...
        inv(eye(num_endog)-ALPHA_ZS'*e_plus), inv(eye(num_exog)-N*e_plus)]);
s_Y_hp{n}=hp(n)^2*s_Y{n};
end



for n = 1 : grid_size
s_yy(:,n)=s_Y{n}(:);
s_yy_hp(:,n)=s_Y_hp{n}(:);
end


s_yy=real(ifft(conj(s_yy')))*2*pi;
s_Y=cell(1,grid_size);


s_yy_hp=real(ifft(conj(s_yy_hp')))*2*pi;
s_Y_hp=cell(1,grid_size);



for n = 1 : grid_size
s_temp = reshape(s_yy(n,:),(1)*num_endog+(1)*num_exog,(1)*num_endog+(1)*num_exog);
s_Y{n}= s_temp;
s_hp_temp = reshape(s_yy_hp(n,:),(1)*num_endog+(1)*num_exog,(1)*num_endog+(1)*num_exog);
s_Y_hp{n}=s_hp_temp;
end


for j=1:length(s_Y)+1
if j-length(s_Y)/2<=0
Gamma_1_spec{j}=s_Y{j+length(s_Y)/2};
Gamma_1_spec_hp{j}=s_Y_hp{j+length(s_Y)/2};
else
Gamma_1_spec{j}=s_Y{j-length(s_Y)/2};
Gamma_1_spec_hp{j}=s_Y_hp{j-length(s_Y)/2};
end
end


if exist('DISPLAY_STD','var')==0 || DISPLAY_STD==1
Standard_Deviations=[char(VARIABLE_NAMES),repmat(char(32),[num_endog+num_exog,3]), num2str(diag(Gamma_1_spec{length(s_Y)/2+1}).^(1/2),'% 0.5f')]
disp(' ')
Standard_Deviations_HP=[char(VARIABLE_NAMES),repmat(char(32),[num_endog+num_exog,3]), num2str(diag(Gamma_1_spec_hp{length(s_Y)/2+1}).^(1/2),'% 0.5f')]
disp(' ')
end

for j=1:length(Gamma_1_spec)
for i=1:num_endog+num_exog
Xcorr_spec{i}(:,j)=Gamma_1_spec{j}(:,i)./((Gamma_1_spec{length(s_Y)/2+1}(i,i).^(1/2))*(diag(Gamma_1_spec{length(s_Y)/2+1}).^(1/2)));
Xcorr_spec_hp{i}(:,j)=Gamma_1_spec_hp{j}(:,i)./((Gamma_1_spec_hp{length(s_Y)/2+1}(i,i).^(1/2))*(diag(Gamma_1_spec_hp{length(s_Y)/2+1}).^(1/2)));
end
end

%%%%%%Begin correlation plots%%%%%%%%%%%
if exist('PLOT_CORRELATION','var')==0 || PLOT_CORRELATION==1 %If either the option has not been set or if the option PLOT_CORRELATION has been set to 1

if exist('CORRELATION_SELECT','var')==0
CORRELATION_SELECT=1:(num_endog+num_exog);
end
if exist('CORRELATION_HORIZON','var')==0
CORRELATION_HORIZON=6;
end

for j=1:length(CORRELATION_SELECT)
    clf('reset')
    set(gcf,'DefaultAxesColorOrder',[hsv(length(CORRELATION_SELECT)-1);0.25, 0.25,0.25])
    CROSS=plot(-CORRELATION_HORIZON:CORRELATION_HORIZON,Xcorr_spec{CORRELATION_SELECT(j)}(CORRELATION_SELECT,grid_size/2+1-CORRELATION_HORIZON:grid_size/2+1+CORRELATION_HORIZON)',':.','MarkerEdgeColor','auto','MarkerSize',8);
    eval(sprintf('title(''Frequency-Domain Cross-Correlations of Selected Variables at t+j with %s at t'')',VARIABLE_NAMES{CORRELATION_SELECT(j)}))
    legend(VARIABLE_NAMES(CORRELATION_SELECT,1),'Location','Best');
    ylabel('Correlation Coefficient');
    xlabel('j');
    axis([-Inf,Inf,-1,1])
    pause;
end
for j=1:length(CORRELATION_SELECT)
    clf('reset')
    set(gcf,'DefaultAxesColorOrder',[hsv(length(CORRELATION_SELECT)-1);0.25, 0.25,0.25])
    CROSS=plot(-CORRELATION_HORIZON:CORRELATION_HORIZON,Xcorr_spec_hp{CORRELATION_SELECT(j)}(CORRELATION_SELECT,grid_size/2+1-CORRELATION_HORIZON:grid_size/2+1+CORRELATION_HORIZON)',':.','MarkerEdgeColor','auto','MarkerSize',8);
    eval(sprintf('title(''HP-Filtered Frequency-Domain Cross-Correlations of Selected Variables at t+j with %s at t'')',VARIABLE_NAMES{CORRELATION_SELECT(j)}))
    legend(VARIABLE_NAMES(CORRELATION_SELECT,1),'Location','Best');
    ylabel('Correlation Coefficient');
    xlabel('j');
    axis([-Inf,Inf,-1,1])
    pause;
end
end
end