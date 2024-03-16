%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%limiting_matrix.m
%
%Used when inifite past expectations are present
%
%This program checks for the existence of user-defined limiting matrix sums
%If they are not present, the program checks for the availability of
%the MATLAB tookbox Symbolic Toolbox and then uses the toolbox to symbolically calculate the limiting
%matrix sums.
%
%If the Symbolic Toolbox is not installed, the user must provide the
%limiting matrix sums him/herself or set the counter it_max_value to some
%finite number to approximate the infinite expectations
%
%THIS VERSION: 1.1 December 9, 2009
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for matj=[65,66,67,70,71]
    if eval(sprintf('eval([''exist(''''%s_inf'''')''])',char(matj)))==0
        eval(sprintf('disp(''%s_inf not provided by user'')',char(matj)))
        if isempty(ver('symbolic'))==0
            eval(sprintf('disp(''Symbolic Toolbox is installed; attempting to create %s_inf symbolically'')',char(matj)))
            eval(sprintf('disp(''If this fails or if the program crashes, please calculate %s_inf manually'')',char(matj)))
            disp(' ')
            disp('Manual calculation should greatly increase performance')
            disp(' ')
            eval(sprintf('eval([''syms %s;''])',num2str(it_name)))
                if eval(sprintf('eval([''isequal(eval(%s_j),zeros(size(eval(%s_j))))==1''])',char(matj),char(matj)))
                    eval(sprintf('eval([''%s_inf=%s_0;''])',char(matj),char(matj)))
                else
                    eval(sprintf('eval([''%s_inf=eval(symsum(eval(%s_j), %s, 1,inf)+%s_0);''])',char(matj),char(matj),num2str(it_name),char(matj)))
                end
        else
            eval(sprintf('disp(''Symbolic Toolbox is not installed, please calculate %s_inf manually'')',char(matj)))
            disp(' ')
        end
    end
end
        

       
 