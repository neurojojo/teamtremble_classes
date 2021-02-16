%Download https://drive.google.com/open?id=1rB26y0EkQ7c-LFinIVVAVvj_Cgd06xcf
%Place ActivePassiveAnalysis folder in 'C:\ActivePassiveAnalysis\'
clear all

addpath( genpath('D:\teamtremble_classes_Lenovo') )
mymice = findDisplacedMice( 'D:\ActivePassiveAnalysis' );

mymice.associateKey('D:\\teamtremble_figures\\cs_all_indices_key.csv');
mymice.loadDisplacements_allMice();

% Manually curated file with three columns (Animal ID, third index of
% displacement for CS+, third index of displacement for CS-)
%
mymice.associateIngressTimes('D:\\teamtremble_figures\\ingress_times_noingress_set_to_10000.csv');

mymice.getIngressScore_allMice();
mymice.getTrembleScore_allMice();

%% Getting an individual mouse analyzed

PW31=displacedMouse('PW31','D:\\ActivePassiveAnalysis\\Other\\PW31\\');
PW31.loadDisplacements();

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

%%

figure; subplot(3,1,1); plot( mymice.mouseObjs.PW31.displacement_csplus(4,:) )
 subplot(3,1,2); plot( mymice.mouseObjs.PW31.fitlines_CSplus(4,:) )
  subplot(3,1,3); plot( mymice.mouseObjs.PW31.corrected_displacement_csplus(4,:) )

  
 
%% CWT NOTES
% Run the cwt manually
% Start off looking at PO372
% Set the debug to the line before the cwt
%
% Run the cwt manually
close all;
figure;
cwt( mymice.mouseObjs.PW17.displacement_csplus( 1 , :), 1000, 'NumOctaves', 10, 'VoicesPerOctave', 48 )
pause(1)
t=get(gca,'Children')
y = t(3).CData;
% Notice that the size of y is: 481 rows
get(gca,'YLim')
% Notice the y-limits of this graph are 0.4239  434.1245
% Take the log values of these, which are the limits of the graph
ylims = log( get(gca,'YLim') )
yrange = linspace( ylims(1), ylims(2), size(y,1) )
[~,closest_to_4hz] = min( abs(exp(yrange)-4) )
[~,closest_to_10hz] = min( abs(exp(yrange)-10) )

[closest_to_4hz,closest_to_10hz] = deal( size(y,1) - closest_to_4hz, size(y,1) - closest_to_10hz )
cwt_ = y( [closest_to_10hz:closest_to_4hz], : );
%cwt_ = y .* conj( y );
%cwt_ = cwt_( [closest_to_10hz:481] , : );
figure; imagesc(cwt_)