%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%matrix_function.m
%
%Used when  past expectations are present
%
%This program uses the strings A_j,B_j,C_j,F_j,G_j, and it_name to
%construct MATLAB anonymous functions for calculations involving
%A_j,B_j,C_j,F_j and G_j for any j
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
eval(sprintf('eval([''target=%s_j;''])',char(matj)))
eval(sprintf('eval([''%s_j=@(%s)%s;''])',char(matj),it_name,target))
end