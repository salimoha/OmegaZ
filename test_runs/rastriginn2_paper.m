function [y,dy] = rastriginn2_paper(x,n)
% rastriginn function 2D example in delta dogs paper
%global n
if nargin <2
    n=2;
end
x=2*(x-0.7);
% n=2;
A=2;
y = A * n * ones(1,size(x,2));
for ii = 1 : 1 : n
    y = y + (x(ii,:) .^ 2 - A * cos(2 * pi * x(ii,:)));
end
% derivative of rasriginn
% dc/dx = 2x+4pi*sin(2pix)
dy = A * n * ones(1,size(x,2));
for ii = 1 : 1 : n
    dy = dy + (2*x(ii,:)  + 2 * pi* A * sin(2 * pi * x(ii,:)));
end
% fun eval
y=y/1.5;
% derivative eval 
dy = dy/1.5;

end
