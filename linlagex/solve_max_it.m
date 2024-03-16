function max_it=solve_max_it(A_j,B_j,C_j,F_j,G_j,tolerance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%solve_max_it.m
%
%Used when infinite past expectations are present
%
%This program uses the anonymous functions A_j,B_j,C_j,F_j and G_j to find
%max_it (I(delta)_max in the paper) such that 
%n>max_it : |(sum_{i=0}^nM_i)_{k,l}-(sum_{i=0}^infinityM_i)_{k,l}|<tol for
%M=A,B,C,F,G;k=1:num_endog;l=1:num_exog
%
%In words: the program finds an iteration after which the matrix sums have
%converged to their, respective, limiting matrices with respect to the
%tolerance criterium 
%
%THIS VERSION: 1.1 December 9, 2009
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


max_it=0;
CON_MAT=@(max_it)max(max(abs([A_j(max_it),B_j(max_it),C_j(max_it),F_j(max_it),G_j(max_it)])));

while CON_MAT(max_it) > tolerance
    max_it=max_it+1;
end
