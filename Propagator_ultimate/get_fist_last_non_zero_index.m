function [first, last]= get_fist_last_non_zero_index(data)
%{
Find the first and last non zero elements in a matrix

by Sergio Bonaque-Gonzalez, PhD. Optical Engineer
sergiob@wooptix.com
July,2019 - Wooptix S.L.
%}

suma=sum(data);

ceros=find(suma~=0);
if ~isempty(ceros)
    first=ceros(1);
    last=ceros(end);
else
    first=1;
    last=length(data);
end
    

    
