function plot_search_const(inter_par_p,inter_par_g,K,xi)
%  keyboard
xv = 0:0.05:1;
tri = delaunayn(xi.');
for ii=1:length(xv)
    for jj=1:length(xv)
        x=[xv(ii); xv(jj)];
%         e=R2-norm(x-xc)^2;
        e = direct_uncer(x,xi,inter_par_p,inter_par_g,tri);
        S(jj,ii)=interpolate_val(x,inter_par_p)-K*e;
        g1(jj,ii)=interpolate_val(x,inter_par_g{1})-K*e;
        g2(jj,ii)=interpolate_val(x,inter_par_g{2})-K*e;
    end
end
figure(3);clf;
contourf(xv,xv,S)
   brighten(0.5)
colormap('bone')
% colorbar
figure(4);clf;
   brighten(0.5)
colormap('bone')
% contourf(xv,xv,max(g1,g2),-1:1:0)
contourf(xv,xv,-max(g1,g2), 0:10:10,  'linestyle', 'none' );
hold on
% plot(xi(1,:),xi(2,:),'ks','MarkerSize',10)
plot(xi(1,:),xi(2,:),'ks',  'MarkerFaceColor','k','MarkerSize', 15) 
end