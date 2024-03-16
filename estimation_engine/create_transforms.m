function [transforms inverse_transforms]=create_transforms(Pparameter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%create_tranforms.m
%
%This program:
%               Creates functions based on the prior distributions to
%               tranform the domains of the parameters to be unrestricted
%Input:         Pparameter: The third and fourth rows correspond to the min
%                    and max of the parameters' domains.
%Ouput:         transforms: a function that transforms the unrestricted domain back to the
%                   target domain
%               inverse_transforms: a function that transforms the target domain into the
%                   unrestricted domain
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

for i=1:size(Pparameter,2)
    if isfinite(Pparameter(3,i))
        if isfinite(Pparameter(4,i))
            transform{1,i}=sprintf('%.20e+(%.20e-%.20e)/(1+exp(x(%d))),',Pparameter(3,i),Pparameter(4,i),Pparameter(3,i),i);
            inverse_transform{1,i}=sprintf('log((%.20e-x(%d))/(x(%d)-%.20e)),',Pparameter(4,i),i,i,Pparameter(3,i));
        else
            transform{1,i}=sprintf('%.20e+exp(x(%d)),',Pparameter(3,i),i);
            inverse_transform{1,i}=sprintf('log(x(%d)-%.20e),',i,Pparameter(3,i));
        end
    else
        if isfinite(Pparameter(4,i))
            transform{1,i}=sprintf('%.20e-exp(-x(%d)),',Pparameter(4,i),i);
            inverse_transform{1,i}=sprintf('-log(%.20e-x(%d)),',Pparameter(4,i),i);
        else
            transform{1,i}=sprintf('x(%d),',i);
            inverse_transform{1,i}=sprintf('x(%d),',i);
        end
    end
end
eval(sprintf('eval([''transforms=@(x)[%s];''])',[transform{:}]))
eval(sprintf('eval([''inverse_transforms=@(x)[%s];''])',[inverse_transform{:}]))