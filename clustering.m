i=10

figure; 

imagesc( mymice.mouseObjs.PO372.cwt_csplus{i}([1:120],:), [0, 0.01] ); hold on; plot( rescale( mean( mymice.mouseObjs.PO372.cwt_csplus{i}([1:120],:), 1 ), 120 ), 'w', 'LineWidth', 2 )
y = mymice.mouseObjs.PO372.cwt_csplus{i}([1:120],:);
y_ = imfilter(y,fspecial('gaussian',5,2) );
figure; hist( log(flat(y_)) )


%%


mytable = summarize(mymice)

output1 = sum(table2array(varfun( @(x) isnan(x), mytable, 'InputVariables', [4:8], 'outputformat', 'table' )),2);
mytable_nonans = mytable( find(output1==0), : );

%

cols_ = [5:7];

mytable_nonans_csplus = mytable_nonans( mytable_nonans.TrialType_Generic=='CS+', : );
mytable_nonans_csplus = sortrows(mytable_nonans_csplus,'Mice');
clustering_table = zscore( table2array(mytable_nonans_csplus(:,cols_)) );

Z = linkage( clustering_table ,'complete');
T = cluster(Z, 'maxclust', 9 );
[~,idx] = sort(T);

lbl = mytable_nonans_csplus.Mice;
%figure; imagesc( clustering_table(idx,:))
lbl = lbl(idx);
%txt_ = arrayfun( @(i) text(1.5,i,lbl{i},'FontSize',6), [1:numel(lbl)] );

figure; 
imagesc( clustering_table(idx,[1,3,2]),[-2,2]);
xticklabels = arrayfun( @(x) regexprep(x,'_',' '), mytable.Properties.VariableNames(cols_([1,3,2])) );
set(gca,'XTick',[1:5],'XTickLabel',xticklabels,'YTick',[1:numel(lbl)],'TickDir','out','YTick',[],'YColor','w'); box off
lbl = lbl(idx);
ylim([0,200])

%% csminus


mytable = summarize(mymice)

output1 = sum(table2array(varfun( @(x) isnan(x), mytable, 'InputVariables', [4:8], 'outputformat', 'table' )),2);
mytable_nonans = mytable( find(output1==0), : );

%

cols_ = [5:7];

mytable_nonans_csplus = mytable_nonans( mytable_nonans.TrialType_Generic=='CS-', : );
mytable_nonans_csplus = sortrows(mytable_nonans_csplus,'Mice');
clustering_table = zscore( table2array(mytable_nonans_csplus(:,cols_)) );

Z = linkage( clustering_table ,'complete');
T = cluster(Z, 'maxclust', 9 );
[~,idx] = sort(T);

lbl = mytable_nonans_csplus.Mice;
%figure; imagesc( clustering_table(idx,:))
lbl = lbl(idx);
%txt_ = arrayfun( @(i) text(1.5,i,lbl{i},'FontSize',6), [1:numel(lbl)] );

figure('color','w'); 
imagesc( clustering_table(idx,[1,3,2]),[-2,2]);
xticklabels = arrayfun( @(x) regexprep(x,'_',' '), mytable.Properties.VariableNames(cols_([1,3,2])) );
set(gca,'XTick',[1:5],'XTickLabel',xticklabels,'YTick',[1:numel(lbl)],'TickDir','out','YTick',[],'YColor','w'); box off
lbl = lbl(idx);
ylim([2,90])

%%


close all
clusters = [7];
figure; imagesc( table2array( mytable_nonans_csplus(find(T==7),cols_) ), [-2,2] ); colormap(jet(200))

myfig = findobj('Type','Figure');
arrayfun( @(x) set( x, 'Position', [520,60,230,740] ), myfig )
myax = findobj('Type','axes');
arrayfun( @(x) set(myax(x),'XTick',[1:5],'XTickLabel',xticklabels,...
    'XTickLabelRotation',90,'YTick',[1:numel(lbl)],'YTickLabel',mytable_nonans_csplus(find(T==clusters(x)),:).Mice), [1:numel(myax)] )

%%

figure; imagesc( table2array( mytable_nonans_csplus(find(T~=7),cols_) ), [-2,2] ); colormap(jet(200))

myfig = findobj('Type','Figure');
arrayfun( @(x) set( x, 'Position', [520,60,230,740] ), myfig )
myax = findobj('Type','axes');
arrayfun( @(x) set(myax(x),'XTick',[1:5],'XTickLabel',xticklabels,...
    'XTickLabelRotation',90), [1:numel(myax)] )


%%

close all
clusters = [1,2,4,7,8,9];
figure; imagesc( table2array( mytable_nonans_csplus(find(T==1),cols_) ), [-2,2] ); colormap(jet(200))
figure; imagesc( table2array( mytable_nonans_csplus(find(T==2),cols_) ), [-2,2] ); colormap(jet(200))
figure; imagesc( table2array( mytable_nonans_csplus(find(T==4),cols_) ), [-2,2] ); colormap(jet(200))
figure; imagesc( table2array( mytable_nonans_csplus(find(T==7),cols_) ), [-2,2] ); colormap(jet(200))
figure; imagesc( table2array( mytable_nonans_csplus(find(T==8),cols_) ), [-2,2] ); colormap(jet(200))
figure; imagesc( table2array( mytable_nonans_csplus(find(T==9),cols_) ), [-2,2] ); colormap(jet(200))

myfig = findobj('Type','Figure');
arrayfun( @(x) set( x, 'Position', [520,60,230,740] ), myfig )
myax = findobj('Type','axes');
arrayfun( @(x) set(myax(x),'XTick',[1:5],'XTickLabel',xticklabels,...
    'XTickLabelRotation',90,'YTick',[1:numel(lbl)],'YTickLabel',mytable_nonans_csplus(find(T==clusters(x)),:).Mice), [1:numel(myax)] )

%%

mytable_ = mytable( mytable.Ingress_time>1000, : );
figure; scatter( mytable_.Ingress_Mag,mytable_.Tremble_Dur,1000,mytable_.TrialType_Generic,'.' ); xlim([0,1]); set(gca,'Colormap',lines(2)); xlabel('Ingress magnitude'); ylabel('Tremble duration');

