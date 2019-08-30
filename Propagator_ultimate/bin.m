function out = bin(in,f)
% f debe ser entero
[m,n] = size(in);

out = zeros(m/f,n/f);

for i = 1:m/f
    for j = 1:n/f
        chop = in((i-1)*f+1:i*f,(j-1)*f+1:j*f);
        out(i,j) = sum(chop(:));
    end
end
end