rng(10)
clearvars x y
figure('color','w'); 
ax = axes('NextPlot','add','TickDir','out');

for i = 1:100
    x(i,:) = [1:100];
    y(i,:) = accumarray([50:60]',ones(1,11)',[100,1]) .* exprnd(.5,100,1);
end

h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','off','Parent',ax); grid off; box off;
colormap((hot(100)))
h.XBinEdges = [1.5:99.5];
h.YBinEdges = [0:.1:max(h.YBinEdges)];
h.EdgeColor = 'none'

plot( mean(y,1), 'LineWidth', 4 )