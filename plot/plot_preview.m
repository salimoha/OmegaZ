function plot_preview(xi,yi,C )
%plot_preview: plots the fucntin evalution and constraits as a fucntion of 
%            iteration
global ms ConBound
h=figure(13);clf;
subplot(ms+1,1,1)
hold on
plot(yi,'-')
plot([0,length(C{1})], [0,0], 'k-')
   ylim([-0.2,0.2])
hold off
ylabel('C0 ')
for ii = 1:ms
    subplot(ms+1,1,ii+1)
    hold on
    plot(C{ii} - ConBound(ii),'-')
    plot([0,length(C{1})], [0,0], 'k-')
    if ii < 3
    ylim([-0.2,0.2]*10)
    else
       ylim([-0.2,0.2])  
    end
    hold off
    ylabel(strcat('C ', num2str(ii)))
end

end

