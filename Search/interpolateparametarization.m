function inter_par= interpolateparametarization(xi,yi,inter_method)
% Calculate the interpolating function

%% polyharmonic spline: inter_method=1. inter_par{1}=inter_method,
% inter_par{2}= w (column vector), inter_par{3}= v (column vector), inter_par{4}=xi.



global bnd1 bnd2



n=size(xi,1);
if inter_method==1
    N = size(xi,2); A = zeros(N,N);
for ii = 1 : 1 : N
    for jj = 1 : 1 : N
        A(ii,jj) = ((xi(:,ii) - xi(:,jj))' * (xi(:,ii) - xi(:,jj))) ^ (3 / 2);
    end
end
V = [ones(1,N); xi];
A = [A V'; V zeros(n+1,n+1)];
% keyboard
% wv = A \ [yi.'; zeros(n+1,1)]; % solve the associated linear system
wv = pinv(A) * [yi.'; zeros(n+1,1)];  % solve the associated linear system
inter_par{1}=1;
inter_par{2} = wv(1:N); inter_par{3} = wv(N+1:N+n+1); 
inter_par{4}= xi;
end

if inter_method==2
[ymin,ind]=min(yi);
xmin=xi(:,ind);
xi(:,ind)=[];
yi(:,ind)=[];
X=xi-repmat(xmin,1,size(xi,2));
Xcl=zeros(n,n); Xul=zeros(n,2);
for ii=1:n
Xcl(ii,ii)=bnd1(ii)-xmin(ii); Xcu(ii,ii)=bnd2(ii)-xmin(ii);
end
A=[X' (X.^2)' zeros(size(X,2),n); -X' -(X.^2)' zeros(size(X,2),n); zeros(n,n) -eye(n) zeros(n,n); -Xcl' -(Xcl.^2)' eye(n); -Xcu' -(Xcu.^2)' eye(n,n)];
B=[yi.'-ymin; 0*yi'; zeros(3*n,1)];
f=[zeros(1,2*n) -ones(1,n)];
options=optimoptions('linprog','Display','none','Algorithm','dual-simplex','MaxIter',100);
x = linprog(f,A,B,[],[],[],[],[],options);
inter_par{1}=2;
inter_par{2}=x(n+1:2*n);
inter_par{3}=x(1:n);
inter_par{4}=xmin;
inter_par{5}=ymin;
end
end