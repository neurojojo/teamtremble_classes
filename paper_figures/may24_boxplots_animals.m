%% 

mytable = summarize( mymice )

%% CS- early versus late trials
mytable = mytable( mytable.TrialType_Generic=='CS-', : );

sortparam = 'Ingress_Mag'
fun_sortparam = sprintf('Fun_%s',sortparam)
figure('color','w');
output1 = varfun(@(x) prctile(x,70), mytable, 'InputVariables', {sortparam}, 'GroupingVariables', {'Mice'} );
[a,b,output1.sorted] = unique( output1.(fun_sortparam) );


tmp_table = innerjoin( mytable, output1 );

% Top plot %
param = 'Ingress_Mag';
lbls = (sortrows(output1,'sorted')); lbls = lbls.Mice;
ax = subplot(3,1,1,'nextplot','add','tickdir','out');
boxplot( tmp_table.(param), tmp_table.sorted )
set(ax,'box','off');
hold on
scatter( tmp_table(tmp_table.Trials<6,:).sorted, tmp_table(tmp_table.Trials<6,:).Ingress_Mag, 'filled', 'jitter', 'on', 'markerfacealpha', 0.5 )
scatter( tmp_table(tmp_table.Trials<6,:).sorted, tmp_table(tmp_table.Trials>=6,:).Ingress_Mag, 'filled', 'jitter', 'on', 'markerfacealpha', 0.5 )
set(gca,'XTickLabel', lbls, 'XTick', [1:numel(lbls)], 'XTickLabelRotation', 90,'ylim',[0,1] )
ylabel( strrep(param,'_','') )

param = 'Tremble_Dur';
% Middle plot %
ax = subplot(3,1,2,'nextplot','add','tickdir','out');
boxplot( tmp_table.(param), tmp_table.sorted )
set(ax,'box','off');
hold on
scatter( tmp_table(tmp_table.Trials<6,:).sorted, tmp_table(tmp_table.Trials<6,:).(param), 'filled', 'jitter', 'on', 'markerfacealpha', 0.5 )
scatter( tmp_table(tmp_table.Trials<6,:).sorted, tmp_table(tmp_table.Trials>=6,:).(param), 'filled', 'jitter', 'on', 'markerfacealpha', 0.5 )
set(gca,'XTickLabel', lbls, 'XTick', [1:numel(lbls)], 'XTickLabelRotation', 90,'ylim',[0,2] )
ylabel( strrep(param,'_','') )


param = 'Tremble_Mag';
% Bottom plot %
ax = subplot(3,1,3,'nextplot','add','tickdir','out');
boxplot( tmp_table.(param), tmp_table.sorted )
set(ax,'box','off');
hold on
scatter( tmp_table(tmp_table.Trials<6,:).sorted, tmp_table(tmp_table.Trials<6,:).(param), 'filled', 'jitter', 'on', 'markerfacealpha', 0.5 )
scatter( tmp_table(tmp_table.Trials<6,:).sorted, tmp_table(tmp_table.Trials>=6,:).(param), 'filled', 'jitter', 'on', 'markerfacealpha', 0.5 )
set(gca,'XTickLabel', lbls, 'XTick', [1:numel(lbls)], 'XTickLabelRotation', 90,'ylim',[0,5] )
ylabel( strrep(param,'_','') )

lines_ = findobj(gcf,'Type','line');
arrayfun( @(x) set(x,'color','k'), lines_ )
delete(findobj('Tag','Outliers'))





%% CS+ versus CS- (all trials)

mytable = summarize( mymice )

sortparam = 'Ingress_Mag'
fun_sortparam = sprintf('Fun_%s',sortparam)
output1 = varfun(@(x) prctile(x,70), mytable, 'InputVariables', {sortparam}, 'GroupingVariables', {'Mice'} );
[a,b,output1.sorted] = unique( output1.(fun_sortparam) );

tmp_table = innerjoin( mytable, output1 );

figure('color','w');
% Top plot %
param = 'Ingress_Mag';
lbls = (sortrows(output1,'sorted')); lbls = lbls.Mice;
ax = subplot(3,1,1,'nextplot','add','tickdir','out');
%boxplot( tmp_table.(param), tmp_table.sorted )
set(ax,'box','off');
hold on
condition = or( tmp_table.TrialType=='CS+1', tmp_table.TrialType=='CS+2' );
scatter( tmp_table(condition,:).sorted, tmp_table(condition,:).Ingress_Mag, 'filled', 'jitter', 'on', 'markerfacealpha', 0.5 );

condition = or( tmp_table.TrialType=='CS-1', tmp_table.TrialType=='CS-2' );
scatter( tmp_table(condition,:).sorted, tmp_table(condition,:).Ingress_Mag, 'filled', 'jitter', 'on', 'markerfacealpha', 0.5 );

set(gca,'XTickLabel', lbls, 'XTick', [1:numel(lbls)], 'XTickLabelRotation', 90,'ylim',[0,1] )
ylabel( strrep(param,'_','') )

param = 'Tremble_Dur';
% Middle plot %
ax = subplot(3,1,2,'nextplot','add','tickdir','out');
%boxplot( tmp_table.(param), tmp_table.sorted )
set(ax,'box','off');
hold on
scatter( tmp_table.sorted, tmp_table.(param), 'filled', 'jitter', 'on', 'markerfacealpha', 0.5 )
scatter( tmp_table.sorted, tmp_table.(param), 'filled', 'jitter', 'on', 'markerfacealpha', 0.5 )
set(gca,'XTickLabel', lbls, 'XTick', [1:numel(lbls)], 'XTickLabelRotation', 90,'ylim',[0,2] )
ylabel( strrep(param,'_','') )


param = 'Tremble_Mag';
% Bottom plot %
ax = subplot(3,1,3,'nextplot','add','tickdir','out');
%boxplot( tmp_table.(param), tmp_table.sorted )
set(ax,'box','off');
hold on
scatter( tmp_table.sorted, tmp_table.(param), 'filled', 'jitter', 'on', 'markerfacealpha', 0.5 )
scatter( tmp_table.sorted, tmp_table.(param), 'filled', 'jitter', 'on', 'markerfacealpha', 0.5 )
set(gca,'XTickLabel', lbls, 'XTick', [1:numel(lbls)], 'XTickLabelRotation', 90,'ylim',[0,5] )
ylabel( strrep(param,'_','') )

lines_ = findobj(gcf,'Type','line');
arrayfun( @(x) set(x,'color','k'), lines_ )
delete(findobj('Tag','Outliers'))





%% Kmeans

output1 = varfun( @(x) nanmean(x), mytable, 'GroupingVariables',{'Mice'},'InputVariables',[5:7],'OutputFormat','table');
output1_array = table2array(output1(:,[3:end]));
out1put1_array_z = zscore(output1_array);
[a,b] = kmeans(  out1put1_array_z,2);
[~,a_sort] = sort( a );
figure; imagesc(  output1_array(a_sort,:) ,[0,1])

%%
figure('color','w'); 
output1_array = table2array( output1(:,[4:6]) );
output1_array_z = zscore( output1_array );
Z = linkage( output1_array_z ,'complete');
dist = pdist( output1_array_z );
leafOrder = optimalleaforder(Z,dist)

%T = cluster(Z, 'maxclust', 4 );
%[~,T_sort] = sort(T);
ax_ = subplot(1,2,1); imagesc(  flipud(output1_array(leafOrder,:)) ,[0,3])
ax__ = subplot(1,2,2); dendrogram(Z,'Reorder',leafOrder); camroll(-90)
set(ax__,'Position',[0.5,.11,.33,.815],'YTick',[1:20])
set(ax_,'YTick',[1:20])

%%


mytable_ = mytable;
pie( mytable_(mytable_.Ingress_Mag>0.84,:).TrialType_Generic )

%%
mytable_ = mytable( mytable.Ingress_onset>0.2, : );
pie( mytable_(or( mytable_.Ingress_Mag>0.83, mytable_.Tremble_Dur>0.83 ),:).TrialType_Generic )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Very good figure for ROC %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mytable = summarize( mymice );
mytable = mytable( mytable.Trials<6, : );

output_csplus = varfun( @(x) nanmean(x), mytable(mytable.TrialType_Generic=='CS+',:), 'InputVariables', {'Ingress_Mag'}, 'GroupingVariables', {'Mice'}, 'OutputFormat', 'table' )
output_csminus = varfun( @(x) nanmean(x), mytable(mytable.TrialType_Generic=='CS-',:), 'InputVariables', {'Ingress_Mag'}, 'GroupingVariables', {'Mice'}, 'OutputFormat', 'table' )

figure('color','w');
subplot(1,2,1); scatter( output_csminus.Fun_Ingress_Mag, output_csplus.Fun_Ingress_Mag, 200, output_csplus.Mice, 'filled' ); title('Ingress')

arrayfun( @(x) text(output_csminus.Fun_Ingress_Mag(x)+0.03, output_csplus.Fun_Ingress_Mag(x), output_csplus.Mice(x) ), [1:20] )

output_csplus = varfun( @(x) nanmean(x), mytable(mytable.TrialType_Generic=='CS+',:), 'InputVariables', {'Tremble_Dur'}, 'GroupingVariables', {'Mice'}, 'OutputFormat', 'table' )
output_csminus = varfun( @(x) nanmean(x), mytable(mytable.TrialType_Generic=='CS-',:), 'InputVariables', {'Tremble_Dur'}, 'GroupingVariables', {'Mice'}, 'OutputFormat', 'table' )
subplot(1,2,2); scatter( output_csminus.Fun_Tremble_Dur, output_csplus.Fun_Tremble_Dur, 200, output_csplus.Mice, 'filled' ); title('Tremble')

arrayfun( @(x) text(output_csminus.Fun_Tremble_Dur(x)+0.03, output_csplus.Fun_Tremble_Dur(x), output_csplus.Mice(x) ), [1:20] )

ax_ = findobj(gcf,'Type','Axes')
arrayfun( @(x) set(x,'XLim',[0,1],'YLim',[0,1]), ax_ )

%% With Tremble Score!

mytable = summarize( mymice );

mytable.Tremble_Score = mytable.Tremble_Mag .* mytable.Tremble_Dur

output_csplus = varfun( @(x) nanmean(x), mytable(mytable.TrialType_Generic=='CS+',:), 'InputVariables', {'Ingress_Mag'}, 'GroupingVariables', {'Mice'}, 'OutputFormat', 'table' )
output_csminus = varfun( @(x) nanmean(x), mytable(mytable.TrialType_Generic=='CS-',:), 'InputVariables', {'Ingress_Mag'}, 'GroupingVariables', {'Mice'}, 'OutputFormat', 'table' )

figure('color','w');
subplot(1,2,1); scatter( output_csminus.Fun_Ingress_Mag, output_csplus.Fun_Ingress_Mag, 200, output_csplus.Mice, 'filled' ); title('Ingress')

arrayfun( @(x) text(output_csminus.Fun_Ingress_Mag(x)+0.03, output_csplus.Fun_Ingress_Mag(x), output_csplus.Mice(x) ), [1:20] )

output_csplus = varfun( @(x) nanmean(x), mytable(mytable.TrialType_Generic=='CS+',:), 'InputVariables', {'Tremble_Score'}, 'GroupingVariables', {'Mice'}, 'OutputFormat', 'table' )
output_csminus = varfun( @(x) nanmean(x), mytable(mytable.TrialType_Generic=='CS-',:), 'InputVariables', {'Tremble_Score'}, 'GroupingVariables', {'Mice'}, 'OutputFormat', 'table' )
subplot(1,2,2); scatter( output_csminus.Fun_Tremble_Score, output_csplus.Fun_Tremble_Score, 200, output_csplus.Mice, 'filled' ); title('Tremble')

arrayfun( @(x) text(output_csminus.Fun_Tremble_Score(x)+0.03, output_csplus.Fun_Tremble_Score(x), output_csplus.Mice(x) ), [1:20] )

ax_ = findobj(gcf,'Type','Axes')
set( ax_(1), 'XLim',[0,2.5],'YLim',[0,2.5],'Nextplot','add' ); plot( ax_(1), [0:2.5],[0:2.5] )
set( ax_(2), 'XLim',[0,1],'YLim',[0,1],'Nextplot','add' ); plot( ax_(2), [0:1],[0:1] )

%% Tremble Ingress correlations


mytable = summarize( mymice );

mytable.Tremble_Score = mytable.Tremble_Dur;

tremble_ingress_cc = arrayfun( @(mouse) mymice.mycorrcoef( mytable( mytable.Mice == mouse , : ).Tremble_Score, mytable( mytable.Mice == mouse , : ).Ingress_Mag ), ...
    unique(mytable.Mice) );


output = varfun( @nanmean, mytable, 'GroupingVariables', {'Mice'}, 'InputVariables', {'Tremble_Score'} );
figure('color','w'); scatter( output.nanmean_Tremble_Score, tremble_ingress_cc, 200, unique(mytable.Mice), 'filled' )
set(gca,'TickDir','out','xlabel',text(0,0,'Tremble Score'),'ylabel',text(0,0,'Tremble-Ingress Correlation') )
xlim([0,0.8])


%% Tremble-ingress correlations by trial

mytable = summarize( mymice );

mytable.Tremble_Score = mytable.Tremble_Mag .* mytable.Tremble_Dur;

mytable_csplus = mytable( mytable.TrialType_Generic=='CS+', : );
tremble_ingress_cc_csplus = arrayfun( @(trial)  mymice.mycorrcoef( mytable_csplus( mytable_csplus.Trials == trial , : ).Tremble_Score, mytable_csplus( mytable_csplus.Trials == trial , : ).Ingress_Mag ), ...
    [1:10] );


mytable_csminus = mytable( mytable.TrialType_Generic=='CS-', : );
tremble_ingress_cc_csminus = arrayfun( @(trial)  mymice.mycorrcoef( mytable_csminus( mytable_csminus.Trials == trial , : ).Tremble_Score, mytable_csminus( mytable_csminus.Trials == trial , : ).Ingress_Mag ), ...
    [1:10] );

figure('color','w'); scatter( [1:10], tremble_ingress_cc_csplus, 200, 'filled' ); hold on;
scatter( [1:10], tremble_ingress_cc_csminus, 200, 'filled' );

%% Tremble-ingress correlations by trial


mytable = summarize( mymice );

mytable.Tremble_Score = mytable.Tremble_Dur;

tremble_ingress_cc = arrayfun( @(trial)  mymice.mycorrcoef( mytable( mytable.Trials == trial , : ).Tremble_Score, mytable( mytable.Trials == trial , : ).Ingress_Mag ), ...
    [1:10] );

figure('color','w'); scatter( [1:10], tremble_ingress_cc, 200, 'filled' ); 

%% Trembly plot

mytable = summarize(mymice);
%figure('color','w'); myscatter = scatter( mytable.Ingress_Fit, mytable.Ingress_Mag, 10*mytable.Trials, mytable.TrialType, 'filled', 'markerfacealpha', 0.8 )
figure('color','w'); 
colors = jet(30); 

ax_ = arrayfun( @(x) subplot(2,4,x),[1:8]);

condition = or( mytable.TrialType=='CS+1', mytable.TrialType=='CS+2' );
myscatter = scatter( ax_(2), log(mytable(condition,:).Tremble_Mag), mytable(condition,:).Tremble_Dur, 50, ones(sum(condition),1), 'filled', 'markerfacealpha', 0.5 )
set(ax_(2),'colormap',colors(30,:),'xlim',[-10,5],'ylim',[0,1], 'tickdir', 'out','xlabel',text(0,0,'Tremble Magnitude'),'ylabel',text(0,0,'Tremble Duration') );

condition = or( mytable.TrialType=='CS-1', mytable.TrialType=='CS-2' );
myscatter = scatter( ax_(4), log(mytable(condition,:).Tremble_Mag), mytable(condition,:).Tremble_Dur, 50, ones(sum(condition),1), 'filled', 'markerfacealpha', 0.5 )
set(ax_(4),'colormap',colors(1,:),'xlim',[-10,5],'ylim',[0,1], 'tickdir', 'out','xlabel',text(0,0,'Tremble Magnitude'),'ylabel',text(0,0,'Tremble Duration') );

for i = 5:8
    pos2 = get(ax_(i),'position');
    set(ax_(i),'position',[pos2(1),0.4,pos2(3),0.1])
    pos2 = get(ax_(i),'outerposition');
end

bins = [-10:1:5];
bins2 = [0:.05:1];
condition = or( mytable.TrialType=='CS+1', mytable.TrialType=='CS+2' );
h(1) = histogram( ax_(1), mytable(condition,:).Tremble_Dur, bins2, 'normalization', 'probability' ); 
set(ax_(1),'YLim',[0,0.4],'View',[-90 90], 'XColor','none','Ycolor','none');
h(2) = histogram( ax_(6), log(mytable(condition,:).Tremble_Mag), bins, 'normalization', 'probability' ); 
set(ax_(6),'YLim',[0,0.25],'YDir','reverse');

condition = or( mytable.TrialType=='CS-1', mytable.TrialType=='CS-2' );
h(3) = histogram( ax_(3), mytable(condition,:).Tremble_Dur, bins2, 'normalization', 'probability' ); 
set(ax_(3),'YLim',[0,0.4],'View',[-90 90], 'XColor','none','Ycolor','none');
h(4) = histogram( ax_(8), log(mytable(condition,:).Tremble_Mag), bins, 'normalization', 'probability' );
set(ax_(8),'YLim',[0,0.25],'YDir','reverse');
arrayfun( @(x) set( x,'box','off'), ax_ )
arrayfun( @(x) set(ax_(x),'xtick',[],'ytick',[], 'xcolor','none','ycolor','none'), [5:8] )
arrayfun( @(x) set(x,'FaceColor',[0.2,0.2,0.2]), h )

set(gcf, 'Position', [160         176        1530         622]);
saveas(gcf, 'may24_scatter_tremble.svg')

% Ingressy plot

mytable = summarize(mymice);
%figure('color','w'); myscatter = scatter( mytable.Ingress_onset, mytable.Ingress_Mag, 10*mytable.Trials, mytable.TrialType, 'filled', 'markerfacealpha', 0.8 )
figure('color','w'); 
colors = jet(30); 

ax_ = arrayfun( @(x) subplot(2,4,x),[1:8]);

condition = and( mytable.Ingress_Mag>=0, or( mytable.TrialType=='CS+1', mytable.TrialType=='CS+2' ));
myscatter = scatter( ax_(2), mytable(condition,:).Ingress_onset, mytable(condition,:).Ingress_Mag, 50, ones(sum(condition),1), 'filled', 'markerfacealpha', 0.5 )
set(ax_(2),'colormap',colors(30,:),'xlim',[0,1], 'ylim', [0,1], 'tickdir', 'out', 'xlabel', text(0,0,'Ingress Onset'), 'ylabel', text(0,0,'Ingress Magnitude') );

condition = and( mytable.Ingress_Mag>=0, or( mytable.TrialType=='CS-1', mytable.TrialType=='CS-2' ));
myscatter = scatter( ax_(4), mytable(condition,:).Ingress_onset, mytable(condition,:).Ingress_Mag, 50, ones(sum(condition),1), 'filled', 'markerfacealpha', 0.5 )
set(ax_(4),'colormap',colors(1,:),'xlim',[0,1], 'ylim', [0,1], 'tickdir', 'out', 'xlabel', text(0,0,'Ingress Onset'), 'ylabel', text(0,0,'Ingress Magnitude') );

for i = 5:8
    pos2 = get(ax_(i),'position');
    set(ax_(i),'position',[pos2(1),0.4,pos2(3),0.1])
    pos2 = get(ax_(i),'outerposition');
end

bins = [0:0.05:1];

condition = and( mytable.Ingress_Mag>=0, or( mytable.TrialType=='CS-1', mytable.TrialType=='CS-2' ));
h(1) = histogram( ax_(3), mytable(condition,:).Ingress_Mag, bins, 'normalization', 'probability' ); 
set(ax_(3),'YLim',[0,0.8],'View',[-90 90], 'XColor','none','Ycolor','none');
h(2) = histogram( ax_(8), mytable(condition,:).Ingress_onset, bins, 'normalization', 'probability' ); 
set(ax_(8),'YLim',[0,0.2],'YDir','reverse');

condition = and( mytable.Ingress_Mag>=0, or( mytable.TrialType=='CS+1', mytable.TrialType=='CS+2' ));
h(3) = histogram( ax_(1), mytable(condition,:).Ingress_Mag, bins, 'normalization', 'probability' ); 
set(ax_(1),'YLim',[0,0.8],'View',[-90 90], 'XColor','none','Ycolor','none');
h(4) = histogram( ax_(6), mytable(condition,:).Ingress_onset, bins, 'normalization', 'probability' );
set(ax_(6),'YLim',[0,0.2],'YDir','reverse');
arrayfun( @(x) set( x,'box','off'), ax_ )
arrayfun( @(x) set(ax_(x),'xtick',[],'ytick',[], 'xcolor','none','ycolor','none'), [5:8] )

arrayfun( @(x) set(x,'FaceColor',[0.2,0.2,0.2]), h )

set(gcf, 'Position', [160         176        1530         622]);
saveas(gcf, 'may24_scatter_ingress.svg')


%% %% Trembly plot (DREADD)

mytable = summarize(mymice_dreadd);
%figure('color','w'); myscatter = scatter( mytable.Ingress_Fit, mytable.Ingress_Mag, 10*mytable.Trials, mytable.TrialType, 'filled', 'markerfacealpha', 0.8 )
figure('color','w'); 
colors = jet(30); 

ax_ = arrayfun( @(x) subplot(2,4,x),[1:8]);

condition = or( mytable.TrialType=='CS+1', mytable.TrialType=='CS+2' );
myscatter = scatter( ax_(2), log(mytable(condition,:).Tremble_Mag), mytable(condition,:).Tremble_Dur, 50, ones(sum(condition),1), 'filled', 'markerfacealpha', 0.5 )
set(ax_(2),'colormap',colors(30,:),'xlim',[-10,5],'ylim',[0,1], 'tickdir', 'out','xlabel',text(0,0,'Tremble Magnitude'),'ylabel',text(0,0,'Tremble Duration') );

condition = or( mytable.TrialType=='CS-1', mytable.TrialType=='CS-2' );
myscatter = scatter( ax_(4), log(mytable(condition,:).Tremble_Mag), mytable(condition,:).Tremble_Dur, 50, ones(sum(condition),1), 'filled', 'markerfacealpha', 0.5 )
set(ax_(4),'colormap',colors(1,:),'xlim',[-10,5],'ylim',[0,1], 'tickdir', 'out','xlabel',text(0,0,'Tremble Magnitude'),'ylabel',text(0,0,'Tremble Duration') );

for i = 5:8
    pos2 = get(ax_(i),'position');
    set(ax_(i),'position',[pos2(1),0.4,pos2(3),0.1])
    pos2 = get(ax_(i),'outerposition');
end

bins = [-10:1:5];
bins2 = [0:.05:1];
condition = or( mytable.TrialType=='CS+1', mytable.TrialType=='CS+2' );
h(1) = histogram( ax_(1), mytable(condition,:).Tremble_Dur, bins2, 'normalization', 'probability' ); 
set(ax_(1),'YLim',[0,0.4],'View',[-90 90], 'XColor','none','Ycolor','none');
h(2) = histogram( ax_(6), log(mytable(condition,:).Tremble_Mag), bins, 'normalization', 'probability' ); 
set(ax_(6),'YLim',[0,0.25],'YDir','reverse');

condition = or( mytable.TrialType=='CS-1', mytable.TrialType=='CS-2' );
h(3) = histogram( ax_(3), mytable(condition,:).Tremble_Dur, bins2, 'normalization', 'probability' ); 
set(ax_(3),'YLim',[0,0.4],'View',[-90 90], 'XColor','none','Ycolor','none');
h(4) = histogram( ax_(8), log(mytable(condition,:).Tremble_Mag), bins, 'normalization', 'probability' );
set(ax_(8),'YLim',[0,0.25],'YDir','reverse');
arrayfun( @(x) set( x,'box','off'), ax_ )
arrayfun( @(x) set(ax_(x),'xtick',[],'ytick',[], 'xcolor','none','ycolor','none'), [5:8] )
arrayfun( @(x) set(x,'FaceColor',[0.2,0.2,0.2]), h )

set(gcf, 'Position', [160         176        1530         622]);
saveas(gcf, 'mymice_dreadd_may24_scatter_tremble.svg')

% Ingressy plot

mytable = summarize(mymice_dreadd);
%figure('color','w'); myscatter = scatter( mytable.Ingress_onset, mytable.Ingress_Mag, 10*mytable.Trials, mytable.TrialType, 'filled', 'markerfacealpha', 0.8 )
figure('color','w'); 
colors = jet(30); 

ax_ = arrayfun( @(x) subplot(2,4,x),[1:8]);

condition = and( mytable.Ingress_Mag>=0, or( mytable.TrialType=='CS+1', mytable.TrialType=='CS+2' ));
myscatter = scatter( ax_(2), mytable(condition,:).Ingress_onset, mytable(condition,:).Ingress_Mag, 50, ones(sum(condition),1), 'filled', 'markerfacealpha', 0.5 )
set(ax_(2),'colormap',colors(30,:),'xlim',[0,1], 'ylim', [0,1], 'tickdir', 'out', 'xlabel', text(0,0,'Ingress Onset'), 'ylabel', text(0,0,'Ingress Magnitude') );

condition = and( mytable.Ingress_Mag>=0, or( mytable.TrialType=='CS-1', mytable.TrialType=='CS-2' ));
myscatter = scatter( ax_(4), mytable(condition,:).Ingress_onset, mytable(condition,:).Ingress_Mag, 50, ones(sum(condition),1), 'filled', 'markerfacealpha', 0.5 )
set(ax_(4),'colormap',colors(1,:),'xlim',[0,1], 'ylim', [0,1], 'tickdir', 'out', 'xlabel', text(0,0,'Ingress Onset'), 'ylabel', text(0,0,'Ingress Magnitude') );

for i = 5:8
    pos2 = get(ax_(i),'position');
    set(ax_(i),'position',[pos2(1),0.4,pos2(3),0.1])
    pos2 = get(ax_(i),'outerposition');
end

bins = [0:0.05:1];

condition = and( mytable.Ingress_Mag>=0, or( mytable.TrialType=='CS-1', mytable.TrialType=='CS-2' ));
h(1) = histogram( ax_(3), mytable(condition,:).Ingress_Mag, bins, 'normalization', 'probability' ); 
set(ax_(3),'YLim',[0,0.8],'View',[-90 90], 'XColor','none','Ycolor','none');
h(2) = histogram( ax_(8), mytable(condition,:).Ingress_onset, bins, 'normalization', 'probability' ); 
set(ax_(8),'YLim',[0,0.2],'YDir','reverse');

condition = and( mytable.Ingress_Mag>=0, or( mytable.TrialType=='CS+1', mytable.TrialType=='CS+2' ));
h(3) = histogram( ax_(1), mytable(condition,:).Ingress_Mag, bins, 'normalization', 'probability' ); 
set(ax_(1),'YLim',[0,0.8],'View',[-90 90], 'XColor','none','Ycolor','none');
h(4) = histogram( ax_(6), mytable(condition,:).Ingress_onset, bins, 'normalization', 'probability' );
set(ax_(6),'YLim',[0,0.2],'YDir','reverse');
arrayfun( @(x) set( x,'box','off'), ax_ )
arrayfun( @(x) set(ax_(x),'xtick',[],'ytick',[], 'xcolor','none','ycolor','none'), [5:8] )

arrayfun( @(x) set(x,'FaceColor',[0.2,0.2,0.2]), h )

set(gcf, 'Position', [160         176        1530         622]);
saveas(gcf, 'mymice_dreadd_may24_scatter_ingress.svg')