function [ALPHA_ZS,BETA_ZS,sorted_eigenvalues,existence_uniqueness]=QZ_solve(A_inf,B_inf,C_inf,F_inf,G_inf,N,num_eqs,num_endog,num_exog,growth_restriction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%QZ_solve.m
%
%Uses Klein (2000) QZ decomposition to find the recursive solution for the
%autonomous recursion of MA-Coefficients. 
%
%THIS VERSION: 1.1 December 9, 2009
%(added section to determine whether imaginary part of solution is
%significant)
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Bring (modified) Binder/Pearson Form into KLEIN FORM
AA=[ zeros(num_eqs,num_endog),-A_inf;eye(num_eqs,num_endog),zeros(num_eqs,num_endog)];
BB=[C_inf,B_inf;zeros(num_eqs,num_endog),eye(num_eqs,num_endog)];
CC=[(F_inf*N+G_inf);zeros(num_eqs,num_exog)];
%QZ Decompostion
[TT, SS, Q,Z] = qz(BB,AA);
%Eigenvalues
all_eigenvalues=ordeig(TT,SS);
%Number of unstable eigenvalues
n_unst=sum(abs(all_eigenvalues)>growth_restriction);
if n_unst==num_endog
    existence_uniqueness=1;
elseif n_unst>num_endog
    existence_uniqueness='unstable';
else
    existence_uniqueness='indeterminate';
end

%Sorts the QZ Decomposition using pre-installed MATLAB implementation of
%LAPACK routines xTGSEN
[TTS,SSS,QS,ZS] = ordqz(TT,SS,Q,Z,abs(all_eigenvalues)<=growth_restriction);
%Sorted eigenvalues
sorted_eigenvalues=diag(TTS)./diag(SSS);

if existence_uniqueness==1
    if rank(ZS(1:2*num_endog-n_unst,1:2*num_endog-n_unst))<num_endog
        existence_uniqueness='non-translatable';
        ALPHA_ZS=[];
        BETA_ZS=[];
    else
        %Assemble final matrix elements
        if num_endog<20 %With a "small" number of endogenous variables, direct calculation quicker
            TEMP_1=kron(N',SSS(2*num_endog-n_unst+1:2*num_endog,2*num_endog-n_unst+1:2*num_endog))...
            -kron(eye(num_exog,num_exog),TTS(2*num_endog-n_unst+1:2*num_endog,2*num_endog-n_unst+1:2*num_endog));
            TEMP_2=QS(2*num_endog-n_unst+1:length(QS),:)*CC;
            vec_Gamma=TEMP_1\reshape(TEMP_2,numel(TEMP_2),1); %Equation 5.7 in Klein (2000)
            Gamma=reshape(vec_Gamma,n_unst,num_exog); %Undo column-vectorization
        else %Otherwise recursive calculation (eqs 5.9-5.13 in Klein (2000)) quicker
            QC=QS(2*num_endog-n_unst+1:length(QS),:)*CC;%define lower half of QS*CC
            for j=1:num_endog
                r_tran=QC(num_endog+1-j,:);%Starting at the end
                    for i=num_endog+2-j:num_endog %and working forwards
                        r_tran=r_tran+TTS(2*num_endog+1-j,i+num_endog)*Gamma(i,:)-SSS(2*num_endog+1-j,i+num_endog)*Gamma(i,:)*N;
                    end
                Gamma(num_endog+1-j,1:num_exog)=r_tran/(SSS(2*num_endog+1-j,2*num_endog+1-j)*N-TTS(2*num_endog+1-j,2*num_endog+1-j)*eye(num_exog));
            end
        end
        %Theta_max_it=ALPHA_ZS*Theta_(max_it-1)+BETA_ZS*(N^(max_it)) (Recursive form for MA coefficients)
        %Using eqs. 5.20 and 5.22 of Klein (2000)
        ALPHA_ZS=ZS(2*num_endog-n_unst+1:2*num_endog,1:2*num_endog-n_unst)/(ZS(1:2*num_endog-n_unst,1:2*num_endog-n_unst));
        BETA_ZS=(ZS(2*num_endog-n_unst+1:2*num_endog,2*num_endog-n_unst+1:2*num_endog)-ZS(2*num_endog-n_unst+1:2*num_endog,1:2*num_endog-n_unst)...
        /(ZS(1:2*num_endog-n_unst,1:2*num_endog-n_unst))*ZS(1:2*num_endog-n_unst,2*num_endog-n_unst+1:2*num_endog))*Gamma;
    end
    else
    ALPHA_ZS=[];
    BETA_ZS=[];
end
if abs(imag([ALPHA_ZS BETA_ZS]))>1E-10
    existence_uniqueness='imaginary';
    ALPHA_ZS=[];
    BETA_ZS=[];
else
    ALPHA_ZS=real(ALPHA_ZS);
    BETA_ZS=real(BETA_ZS);
end