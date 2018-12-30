clear all
close all

for k=1:20
n=20*k;
h=2/(n-1);
y=-1:2/(n-1):1;
T=-2*diag(ones(n,1))+ diag(ones(n-1,1),1)+diag(ones(n-1,1),-1);
T=T/h^2;
U=[zeros(n,n) diag(y); diag(y) zeros(n,n)];
T=[T zeros(n,n); zeros(n,n) T];
T1=T*inv(U);
v=eig(T1);
Re(k)=min(v(v>0))
max(eig(T+Re(k)*U))
end
plot(20:20:400,real(Re))
xlabel('Discreitization size: n','fontsize',20,'fontweight','bold')
ylabel('Re_{Critical}','fontsize',20,'fontweight','bold')
grid on
set(gca,'FontSize',20)


