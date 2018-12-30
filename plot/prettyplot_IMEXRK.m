function  prettyplot_IMEXRK( xi, yi, C )
% plots the constraints, max violation, best point, all ci, and obj. vs iteration 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
global FUN_MAX  CON_MAX
h=figure(13); clf;
set(h,'units','normalized','outerposition',[0 0 0.9 0.9])
%default values
op.fontsize = 12;
op.quiet = 0;
% op.dirpath = conf.outputDir;
op.dirpath = './';
op.width = 1500;
op.height = 900;
op.visible = 1;
op.title = 1;
op.cache = 0;
op.output = 0;
op.embed = 0;
op.save = 0;
op.label = 1; %Settign to zero removes labels on figures
ah(1) = axes('Position', [0.015 0.48 0.3 0.47]);
plot(C{1},'k-','linewidth',3)
hold on
% plot(C{2},'--','linewidth',3)
plot([0,length(C{1})] , [0,0], 'r-.','linewidth',3)

%
hold off
title('C 1-2')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ah(2) = axes('Position', [0.34 0.48 0.30 0.47]);
hold on
plot(yi,'k','linewidth',3)
%
plot(C{3},'-','linewidth',3)
% plot(C{4},'--','linewidth',3)
%
plot([0,length(C{1})] , [0,0], 'r-.','linewidth',3)
hold off
grid on
title('objective function & C 3-4')
% axis off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ah(3) = axes('Position', [0.675 0.48 0.30 0.47]);
% plot(yi,'k','linewidth',3)
hold on
plot(C{5},'k','linewidth',3)
plot(C{6},'linewidth',3)
% plot(C{7},'linewidth',3)
plot([0,length(C{1})] , [0,0], 'r-.','linewidth',3)
ylim([-CON_MAX,CON_MAX])
grid on
title(' radicals and bounding c34')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% second row of plots
ah(4) = axes('Position', [0.015 0.05 0.30 0.3]);
for i = 1:length(yi)
%     [VM0(i),indm(i)] = max([yi(i),C{1}(i),C{2}(i), C{3}(i),C{4}(i),C{5}(i), C{6}(i),C{7}(i) ]);
    [VM0(i),indm(i)] = max([yi(i),C{1}(i),C{2}(i), C{3}(i),C{4}(i),C{5}(i), C{6}(i) ]);
    % VM(i) = max([yi(i),C{1}(i), C{3}(i) ]);
%     [Vm0(i),ind(i)] = min([yi(i),C{1}(i),C{2}(i), C{3}(i),C{4}(i),C{5}(i), C{6}(i),C{7}(i) ]);
    [Vm0(i),ind(i)] = min([yi(i),C{1}(i),C{2}(i), C{3}(i),C{4}(i),C{5}(i), C{6}(i) ]);
     VM(i) = max([yi(i),C{1}(i), C{3}(i) ]);
end
hold on
plot(VM0,'k','linewidth',3)
plot(Vm0,'k--','linewidth',3)
% plot(VM,'b--','linewidth',3)
hold off
grid on
% ylim([0,2])
title('Maximum violation')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ah(5) = axes('Position', [0.35 0.05 0.3 0.3]);
hold on
% trajectory of maximum violation
for ii = 1:length(yi)
    [ym(ii),indx(ii)] = min(VM0(1:ii));
    Xvm(:,ii) = xi(:,indx(ii));
end
hold on
plot(Xvm(1,:),'k-','linewidth',3)
plot(Xvm(2,:),'-','linewidth',3)
plot(Xvm(3,:),'-','linewidth',3)
% plot(VM,'k--','linewidth',3)
hold off
grid on
ylim([0,1])
title('best point trajectory')
set(ah, 'FontSize', op.fontsize);
% positions of xi as iterations
ah(6) = axes('Position', [0.675 0.05 0.3 0.3]);
hold on; grid on
plot(xi(1,:),'k-','linewidth',3)
plot(xi(2,:),'-','linewidth',3)
plot(xi(3,:),'-','linewidth',3)
ylim([0,1])
title('position of c_i')
set(ah, 'FontSize', op.fontsize);


end

