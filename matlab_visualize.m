%% setup graph

% import csv files SiouxFalls_net and SiouxFalls_bis_node manually in
% MATLAB

s = table2array(SiouxFallsnet(:,2));
t = table2array(SiouxFallsnet(:,3));
xdata = table2array(SiouxFallsbisnode(:,2));
ydata = table2array(SiouxFallsbisnode(:,3));

%% import flows
%f = table2array(SiouxFallsbisresult(:,2));
%linktimes = table2array(SiouxFallsbisresult(:,3));

G = cell(1);
numgraphs = size(f,2);
for i = 1:numgraphs
    G{i} = digraph(s,t,f(:,i));
end

%%
allf = reshape(f,[size(f,1)*size(f,2),1]);
fmax = max(allf);

nbins = 100;
colors = parula(nbins);

[~,~,bin] = histcounts(allf,nbins);
fcolors = colors(bin,:);

startind = 1;
nedges = length(f);
for i=1:size(f,2)
    thisf = f(:,i);
    fwidths = max(0.01,10*thisf/fmax);
    h = figure;
    plot(G{i},'Xdata',xdata,'Ydata',ydata,'EdgeColor', fcolors(startind:startind+nedges-1,:), 'EdgeLabel',round(thisf),'LineWidth',fwidths,'EdgeFontsize',10)
    startind = startind+nedges;
    xlabel('X Location','Fontsize',14);
    ylabel('Y Location','Fontsize',14);
    title('Traffic Network Visualization: Sioux Falls with Capacity Constraints','Fontsize',16);
    c = colorbar;
    set(get(c,'Label'),'string','Link flow');
    set(gca,'Fontsize',12);
    
    set(gcf,'Position',[1, 1, 933, 722]);
    set(gca,'Position',[0.1239, 0.1100, 0.7385, 0.8150]);
    annotation('Textbox',[.05 .06 .01 .02],'String',['Link capacity: ',num2str(csol(i))],'FitBoxToText','on','Fontsize',14);
    
    frame(i) = getframe(gcf); 
    
    % for gif creation
%     im = frame2im(frame); 
%     [imind,cm] = rgb2ind(im,256); 
%     
%     % Write to the GIF File 
%     if i == 1 
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%     else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%     end 
end
close 'all'

v = VideoWriter('SiouxFallsnet');
v.FrameRate = 1;
open(v);
writeVideo(v,frame);
close(v);