% Representation of subinterpolation 
clear all
close all
clc
warning('off','all')


c =-.2;
w=10;
fun=@(x) 1*(x-c).^2+sin(w*(x-c));
gfun=@(x) 1*(x-c).^2+sin(w*((x-c)+0.15));
% xi=[0.001 1 0.2 0.3 0.5 0.6 0.12 0.07 .17 .15];
xi=[0.1 0.7 0.3];
xiu=[0 1];
for ii=1:length(xi)
    yi(ii)=fun(xi(:,ii));
    gi(ii)=gfun(xi(:,ii));
end
%T=ones(1,3);
%K=2; L=1;
%[xi,ind]=sort(xi); yi=yi(ind);
xx=0:0.01:1;
for ii=1:length(xx)
yy(ii)=fun(xx(ii));
gg(ii)=gfun(xx(ii));
end
figure(1)
hold on
plot(xx,yy,'k-',xx,gg,'k--','linewidth',2.5)
plot(xx,max(yy,gg)+1.,'k-.','linewidth',3)
legend('f(x)','g(x)','max\{f,g\}')

%%
xiT=[xi xiu]; xiT=sort(xiT);
xc=(xiT(1:end-1)+xiT(2:end))/2;
d=(xiT(2:end)-xiT(1:end-1))/2;
%[inter_par,yp]=interpolateparametarization(xi,yi,sigma0./sqrt(T),1);
[inter_par]=interpolateparametarization(xi,yi,1);
[inter_par_g]=interpolateparametarization(xi,gi,1);
%x=0.01; e=max(d.^2-(x-xc).^2);


%
for ii=1:length(xx)
    x=xx(ii);
    %fr(ii)= funr(x);
    f(ii)=fun(x);
    p(ii) = interpolate_val(x,inter_par);
    pg(ii) = interpolate_val(x,inter_par_g);
    e(ii) = max(d.^2-(x-xc).^2);
    sc(ii)= (max(p(ii),pg(ii))+0.5455)/e(ii);
    %sc(ii) = max(p(ii)-;
    %sn(ii)= subinterpolate_search(x,e(ii),xi,yi,K);
    %sn(ii)=subquad_search(x,e(ii),xi,yi,K);
end
for ii=1:2
    x=xiu(ii);
    pU(ii) = interpolate_val(x,inter_par);
    pgU(ii) = interpolate_val(x,inter_par_g);
    eU(ii) = mindis(x,xi);
     sd(ii)= (max(pU(ii),pgU(ii))+0.5455)/eU(ii);
end


[tc,indc]=min(sc);
x=xx(indc);
pm = interpolate_val(x,inter_par);
pgm = interpolate_val(x,inter_par_g);
em = mindis(x,xi);
sdm= (max(pm,pgm)+0.5455)/em;
     
% %%
figure(2)
plot(xx,pg,'k-',xx,p,'k--', xx, max(p,pg)+2, 'k-.') 
hold on
grid on
plot(xi,yi,'ks', 'Markersize',10, 'MarkerFacecolor', 'k')
 plot(xi,gi,'bp', 'Markersize',10, 'MarkerFacecolor', 'k')
 
figure(3)
plot(xx,1000*e,'r-',xx, sc,'k-') 
hold on
grid on
ylim([0 100])

figure(4)
plot(xiu,sd,'rs',xx(indc), sdm,'ks') 
hold on
grid on
%ylim([0 100])
% % errorbar(xi,yi,sigma0./sqrt(T),'k.','linewidth',2)
% %    sd = min(yp,2*yi-yp)-L*sigma0./sqrt(T);
% %[td,indd]=min(sd);
% %hold on
 
% %%
% %text(xx(indc),tc-.05,'x_k','fontsize',20)
% % plot(xx(indc),tc,'k*','Markersize',10)
% %text(xi(indd),td-.05,'z_k','fontsize',20)
% %set(gca,'YTickLabel',[])
% %set(gca,'XTickLabel',[])
% % ylim([-.3 10])
% grid on
%     
% % 
% % 
% % [t,ind]=min(sn);
% % ind
% % xx(ind)
% % x=xx(ind);
% % ii=ind;
% % [y,abc]=subquad_search(x,e(ii),xi,yi,K)
% % y_inter=abc(1)*xx.^2+ abc(2)*xx+abc(3);
% % hold on
% % plot(xx,y_inter,'r-')