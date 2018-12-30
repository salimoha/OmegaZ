function plot2(x,y)
if nargin <2
plot(x,'-','linewidth',3)
else    
plot(x,y,'-','linewidth',3)
end
set(gca, 'XGrid', 'on', 'YGrid', 'on',  'FontSize', 24);
% set(gcf, 'Position' , [60 60 700 500], 'color' , [1 1 1], 'PaperPositionMode' , 'auto', 'InverthardCopy' , 'off' );
% set(gcf, 'color' , [1 1 1], 'PaperPositionMode' , 'auto', 'InverthardCopy' , 'on' );
set(gcf, 'color' , [1 1 1], 'PaperPositionMode' , 'auto', 'InverthardCopy' , 'on' );
end