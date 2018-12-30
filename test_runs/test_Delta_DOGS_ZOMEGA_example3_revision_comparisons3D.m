% Search.Method=1 constant_K
% Search.Method=2 Adaptive_K

clear all; close all; clc
% global n m ms bnd1 bnd2 Ain bin acon Search r xi tri
global n m ms bnd1 bnd2  Search  Ain bin tri MESH_SIZE iter_max
% Search.constant=0;
%
n=3;
ms=1;
RES = 1e-2;
x_star=ones(n,1) *0.153;

%
Method = 2;
%
x0=[0;0;0]; KCf = 1; KC2 = 1;
fun=@(x) ( (x(1,:)-x0(1)).^2+(x(2,:)-x0(2)).^2 +(x(3,:)-x0(3)).^2 )*KCf;

con{1}=@(x) (rastriginn2(x)-0.5 )*KC2; %the same
% Var = 0.0721;
Var=fun(x_star)+1e-3;

%
MESH_SIZE=8; % grid size% % Mss=20;
% interpolaion strategy
inter_method=1;

% Calculate the initial points
xE=ones(n,1)*0.5;
delta0=0.15;
for ii=1:n
    e=zeros(n,1); e(ii)=1;
    xE(:,ii+1)=xE(:,1)+delta0*e;
end
% calculates acon , yi, C
lob=zeros(n,1); upb=ones(n,1);
Search.method = Method;
Search.constant = Var;
bnd1 = lob;
bnd2 = upb;
xU=bounds(bnd1,bnd2, n);
% Input the equality constraints
Ain=[eye(n);-eye(n)];
bin=[bnd2 ;-bnd1];
% Calculate the function evaluation at initial points
for ii=1:ms
    acon{ii}=[];
end
for ii=1:size(xE,2)
    ind=1:2*n;
    ind=ind(Ain*xE(:,ii)-bin>-0.01);
    for jj=1:length(ind)
        acon{ind(jj)}=[acon{ind(jj)} ii];
    end
    yE(ii)=fun(xE(:,ii));
    for jj=1:ms
        C{jj}(ii)=con{jj}(xE(:,ii));
    end
end
x_prev = xE(:,1);
y_prev = yE(:,1);
delta_tol=0.2;

%
iter_max = 100;
y0=Search.constant;

for kkk=1:16
    for k=1:iter_max
        %             save surr_pts_2 yi xi C tri acon
        disp( [' started iteration ' ,num2str(k), ' ... ', 'with number func. eval = ',num2str(size(xE,2))] )
        tri=delaunayn([xE xU].');
         xi=[xE xU];
        
        %%%%%%%%%%% Modify the interpolations %%%%%%%%%%%%%%%%%%%
        inter_par_p= interpolateparametarization(xE,yE,inter_method);
        for jj=1:ms
            inter_par_g{jj}= interpolateparametarization(xE,C{jj},inter_method);
        end
        
        %%%%%%%%% Calculate Discrete search function %%%%%%%%%%%
        yup=zeros(1,size(xU,2));
        for ii=1:size(xU,2)
            yup(ii)=estimate_max_cons_val(xU(:,ii),inter_par_p,inter_par_g,y0,ms)/mindis(xU(:,ii),xE);
%             yup(ii)=inf;
        end
        
        % Perform the search
        while 1
%             keyboard
            [xm ym(k) Cs(k) indm2]= ...
                tringulation_search_constraints(inter_par_p,inter_par_g,[xE xU], tri)
            %keyboard
            xm=round((xm-bnd1)./(bnd2-bnd1).*MESH_SIZE)./MESH_SIZE.*(bnd2-bnd1)+bnd1;
            
            [xm,xE,xU,newadd,success]=...
                points_neighbers_find(xm,xE,xU);
            
            if success==1
                break
            else
                yup=[yup estimate_max_cons_val(xm,inter_par_p,inter_par_g,y0,ms)/mindis(xm,xE)];
            end
        end
        
        % stopping criteria
        if mindis(xm,[xE xU])<1e-6
            break
        end
        
        if (estimate_max_cons_val(xm,inter_par_p,inter_par_g,y0,ms)/mindis(xm,xE)>min(yup) || mindis(xm,xU)<1e-6)
            [t,ind]=min(yup);
            xm=xU(:,ind); xU(:,ind)=[];
        end
%         % feasible constraint projection
%         if mindis(xm,[xE xU])<1e-6
%             break
%         end
       
        
        
        % Perform the function evalutions
        %     Evaluated set
        xE=[xE xm];
        yE = [yE fun(xm)];
        for jj=1:ms
            C{jj}=[C{jj} con{jj}(xm)];
        end
        
       %
        if Search.method ==2
%             figure(11);
%             subplot(2,1,1)
%             plot(yE)
%             subplot(2,1,2)
%             plot(C{1})
            
                        figure(11);
            subplot(2,1,1)
            plot(yE-Var)
            subplot(2,1,2)
            plot(max(yE-Var, C{1}))
%        
%             figure(12);clf;
%             plot(xE(1,:),xE(2,:),'x')
%             hold on
%             plot(xU(1,:),xU(2,:),'o')
%             plot(xm(1,:),xm(2,:),'ks', 'MarkerSize',10)
%             tt=0:0.01:1;
% %             plot( rastriginn(tt),tt, 'r', 'linewidth', 3)
%             plot(tt, rastriginn(tt), 'r', 'linewidth', 3)
%            
%                 triplot(tri,xi(1,:),xi(2,:)) 
            drawnow
        end
        
        
        
        iplot_check=0;
        if iplot_check==1
            clear U_f p_f G_c C_l Err TT SS G1 G2 xi
            xi=[xE xU];
            xv=0:0.01:1;
            for ii=1:length(xv)
                for jj=1:length(xv)
                    %         U_f(ii,jj)=fun([xv(ii) ;xv(jj)]);
                    P_f(ii,jj)=interpolate_val([xv(ii) ;xv(jj)],inter_par_p);
                    C_l(ii,jj)=con{1}([xv(ii) ;xv(jj)]);
                    %         G_c(ii,jj)=interpolate_val([xv(ii) ;xv(jj)],inter_par_g{1});
                    G1=interpolate_val([xv(ii) ;xv(jj)],inter_par_g{1});
                    G2=interpolate_val([xv(ii) ;xv(jj)],inter_par_g{2});
                    [Err(ii,jj),TT(ii,jj),SS(ii,jj)] = direct_uncer([xv(ii) ;xv(jj)],xi,inter_par_p,inter_par_g,tri);
                    CC(ii,jj) = max([G1, G2])- Search.constant * Err(ii,jj);
                    %          CC_2(ii,jj) = max([G1, G2])- Search.constant * Err(ii,jj);
                    % %        if CC(ii,jj) >0
                    % % %            CC(ii,jj) = nan;
                    % % CC(ii,jj) = 0;
                    % %        else
                    % %            CC(ii,jj) = 1;
                    % %        end
                    
                    
                end
            end
            %
            % h=figure; clf;
            % ah(1) = axes('Position', [0.015 0.48 0.3 0.47]);
            figure(10); clf;
            axis square
            box on
            hold on
            contourf(xv,xv,-CC.', 0:10:10,  'linestyle', 'none' );
            grid on
            grid minor
            plot(x_star(1),x_star(2),'kp','MarkerFaceColor','b', 'MarkerSize', 18)
            brighten(0.5)
            colormap('bone')
            % colorbar
            % caxis([-.010,0])
            % contourf(xv,xv,(CC).')
            % contourf(xv,xv,(CC).','linestyle', 'none',0.1:0.1)
            % set(h,'linestyle','none');
            % caxis([-0.0001,0.0001])
            % colormap(map)
            %
            %   plot(xi(1,:),xi(2,:),'ro', 'MarkerSize', 13)
            if Ex ==1
                plot(tt, b*sin(pi*tt), 'k', 'linewidth', 3)
            elseif Ex==2
                plot(tt, rastriginn(tt), 'k', 'linewidth', 3)
            end
            %% plot the evaluated points
            plot(xi(1,1:4),xi(2,1:4),'ks',  'MarkerFaceColor','k','MarkerSize', 15)
            plot(xi(1,5:end),xi(2,5:end),'ks','MarkerFaceColor','k', 'MarkerSize', 15)
            %   plot(xi(1,4:end-1),xi(2,4:end-1),'ro', 'MarkerSize', 13)
            plot(xi(1,end),xi(2,end),'wx', 'MarkerSize', 15)            
            % caxis([-0.5,0])
            % colorbar
            % title(['iter = ', num2str(k)], 'FontSize', 20)
            % title('max (g_l(x) - K e(x))')
            % axis off
            set(gca, 'ytick', [])
            set(gca, 'xtick', []) 
            drawnow
            
        end        
        %%
        disp([num2str(k/iter_max*100), ' % Completed'])
    end
%     keyboard
    MESH_SIZE=MESH_SIZE*2;
end

%%

save results_OmegaZ_paper3D.mat