function y = rastriginn(x)
% rastriginn function 2D example in delta dogs paper
%global n
x=3.5*(x-0.7);
n=1;
A=2;
y = A * n * ones(1,size(x,2));
for ii = 1 : 1 : n
    y = y + (x(ii,:) .^ 2 - A * cos(2 * pi * x(ii,:)));
end
y=y/12+0.1;
end
