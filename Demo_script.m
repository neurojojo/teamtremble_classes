%Download https://drive.google.com/open?id=1rB26y0EkQ7c-LFinIVVAVvj_Cgd06xcf
%Place ActivePassiveAnalysis folder in 'C:\ActivePassiveAnalysis\'
clear all

addpath( genpath('C:\teamtremble_classes') )
mymice = findDisplacedMice( 'C:\ActivePassiveAnalysis2' );

%
mymice.associateKey('C:\\teamtremble_classes\\cs_all_indices_key.csv');

mymice.loadDisplacements();

mymice.makeList()
% Manually curated file with three columns (Animal ID, third index of
% displacement for CS+, third index of displacement for CS-)


mymice.associateIngressTimes('C:\\teamtremble_classes\\ingress_times_noingress_set_to_10000.csv');

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


