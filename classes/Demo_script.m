%Download https://drive.google.com/open?id=1rB26y0EkQ7c-LFinIVVAVvj_Cgd06xcf
%Place ActivePassiveAnalysis folder in 'C:\ActivePassiveAnalysis\'
clear all

addpath( genpath('C:\teamtremble_classes') )
mymice = findDisplacedMice( 'C:\ActivePassiveAnalysis2' );

%%
mymice.associateKey('C:\\teamtremble_figures\\cs_all_indices_key.csv');

mymice.loadDisplacements_allMice();

% Manually curated file with three columns (Animal ID, third index of
% displacement for CS+, third index of displacement for CS-)
%%
mymice.associateIngressTimes('C:\\teamtremble_figures\\ingress_times_noingress_set_to_10000.csv');

%%
mymice.fitLines_allMice('coeffs');
mymice.fitLines_allMice('fitlines');
mymice.getIngressScore_allMice();
mymice.getTrembleScore_allMice();

%%

csplus_ = mytable( mytable.TrialType_Generic=='CS+',:)
csminus_ = mytable( mytable.TrialType_Generic=='CS-',:)

figure; ax_ = arrayfun( @(x) subplot(10,1,x,'NextPlot','add','XLim',[0,.1],'YLim',[0,.1]) , [1:10] );
colors = gray(10);
arrayfun( @(x) loglog(ax_(x), csplus_( csplus_.Trials==x, : ).Tremble_Mag, csminus_( csminus_.Trials==x, : ).Tremble_Mag, '.', 'color', colors(x,:), 'MarkerSize', 24 ), [1:10] )

%%


figure; ax_ = axes('NextPlot','add', 'XLim', [-4,4], 'YLim', [-4,4], 'XScale', 'log', 'YScale', 'log')
colors = copper(10);
arrayfun( @(x) loglog( csplus_( csplus_.Trials==x, : ).Tremble_Mag, csminus_( csminus_.Trials==x, : ).Tremble_Mag, '.', 'color', colors(x,:), 'MarkerSize', 24 ), [1:10] )

%%


figure; ax_ = axes('NextPlot','add')%,'XScale', 'log', 'YScale', 'log')
colors = copper(10);
arrayfun( @(x) scatter( csplus_( csplus_.Trials==x, : ).Tremble_Dur, csminus_( csminus_.Trials==x, : ).Tremble_Dur, 'ko', 'MarkerFaceColor', colors(x,:) ), [1:10] )

%%

figure; ax_ = axes('NextPlot','add')%,'XScale', 'log', 'YScale', 'log')
colors = copper(10);
arrayfun( @(x) scatter( csplus_( csplus_.Trials==x, : ).Ingress_Mag, csminus_( csminus_.Trials==x, : ).Ingress_Mag, 'ko', 'MarkerFaceColor', colors(x,:) ), [1:10] )

%%


figure; ax_ = axes('NextPlot','add')%,'XScale', 'log', 'YScale', 'log')
colors = copper(10);
arrayfun( @(x) scatter( csplus_( csplus_.Trials==x, : ).Ingress_Fit, csminus_( csminus_.Trials==x, : ).Ingress_Fit, 'ko', 'MarkerFaceColor', colors(x,:) ), [1:10] )


%%

result=mymice.search('m[12]m[12]m[12]')
mymice.mouseObjs.PW17.trembleScore_csplus.Dur


%% Compare end of tremble to beginning of onset

% Four heatmaps are produced: 
% mean_image_csplus
% mean_image_csplus_above_thr
% mean_image_csminus
% mean_image_csminus_above_thr

% Mean of CWTs for CS+
cwt_matrix_struct = {};

for i = 1:numel(myfields);
cwt_matrix_struct{i} = nanmean(mymice.mouseObjs.(myfields{i}).ingressAlignedCWT_csplus(),3);
end;

mean_image_csplus = zeros(201,14001);
for i = 1:numel(myfields);
    mean_image_csplus = nansum( cat(3, mean_image_csplus,cwt_matrix_struct{i}) , 3);
end
mean_image_csplus = mean_image_csplus/numel(myfields);



% Number of trials exceeding tremble threshold for CS+
for i = 1:numel(myfields);
cwt_matrix_struct{i} = nansum(mymice.mouseObjs.(myfields{i}).ingressAlignedCWT_csplus()>0.005,3);
end;

mean_image_csplus_above_thr = zeros(201,14001);
for i = 1:numel(myfields);
    mean_image_csplus_above_thr = mean_image_csplus_above_thr + cwt_matrix_struct{i};
end

% Mean of CWTs for CS-
cwt_matrix_struct = {};

for i = 1:numel(myfields);
cwt_matrix_struct{i} = nanmean(mymice.mouseObjs.(myfields{i}).ingressAlignedCWT_csminus(),3);
end;

mean_image_csminus = zeros(201,14001);
for i = 1:numel(myfields);
    mean_image_csminus = nansum( cat(3, mean_image_csminus,cwt_matrix_struct{i}) , 3);;
end
mean_image_csminus = mean_image_csminus/numel(myfields);

% Number of trials exceeding tremble threshold for CS-
for i = 1:numel(myfields);
cwt_matrix_struct{i} = nansum(mymice.mouseObjs.(myfields{i}).ingressAlignedCWT_csminus()>0.005,3);
end;

mean_image_csminus_above_thr = zeros(201,14001);
for i = 1:numel(myfields);
    mean_image_csminus_above_thr = mean_image_csminus_above_thr + cwt_matrix_struct{i};
end
%% Plotting the variables created in the cell above
figure('position',[567   262   822   357],'color','w'); 
ax = arrayfun( @(x) subplot(1,2,x), [1:2] )
imagesc(ax(1), mean_image_csplus, [0,0.005]); colormap( jet(100) );  
imagesc(ax(2), mean_image_csminus, [0,0.005]);
arrayfun( @(x) set(x,'xlabel',text(0,0,'Time before ingress (s)'),'ytick',[80,130],'yticklabel',{'10Hz','4Hz'},'xtick',[0:2000:14000],'xticklabel',[-14:2:0],'xlim',[4000,14000]), ax );
set( ax(1), 'Title', text(0,0,sprintf('CS+ (%i ingresses)',Ntrials_CSplus),'FontWeight','normal','FontSize',14) )
set( ax(2), 'Title', text(0,0,sprintf('CS- (%i ingresses)',Ntrials_CSminus),'FontWeight','normal','FontSize',14) )
t=suptitle('Averaged cwt before ingress')
set(t,'FontSize',18)

figure('position',[567   262   822   357],'color','w'); 
ax = arrayfun( @(x) subplot(1,2,x), [1:2] )
imagesc(ax(1), mean_image_csplus_above_thr, [0,40]); colormap( jet(100) ); 
imagesc(ax(2), mean_image_csminus_above_thr, [0,40]); 
arrayfun( @(x) set(x,'xlabel',text(0,0,'Time before ingress (s)'),'ytick',[80,130],'yticklabel',{'10Hz','4Hz'},'xtick',[0:2000:14000],'xticklabel',[-14:2:0],'xlim',[4000,14000]), ax );
[Ntrials_CSplus,Ntrials_CSminus] = deal( max(max(mean_image_csplus_above_thr)), max(max(mean_image_csminus_above_thr)) );
set( ax(1), 'Title', text(0,0,sprintf('CS+ (%i ingresses)',Ntrials_CSplus),'FontWeight','normal','FontSize',14) )
set( ax(2), 'Title', text(0,0,sprintf('CS- (%i ingresses)',Ntrials_CSminus),'FontWeight','normal','FontSize',14) )
t=suptitle('Number of ingress trials where cwt>threshold for tremble')
set(t,'FontSize',18)
%imagesc(ax(3), mean_image_csplus_above_thr)
%imagesc(ax(4), mean_image_csminus_above_thr)