figure('color','w');
mytable = summarize( mymice );
subplot(1,2,1); scatter( log(mytable.Tremble_Mag(mytable.Ingress_Mag>0.5)), mytable.Tremble_Dur(mytable.Ingress_Mag>0.5), 100, 'filled', 'markerfacealpha',0.5 ); set(gca,'Colormap',[0,0,0;1,0,0]); xlim([-10,5]); xlabel('TrembleMag');ylabel('TrembleDur'); title('Ingress'); ylim([0,1.1])
subplot(1,2,2); scatter( log(mytable.Tremble_Mag(mytable.Ingress_Mag<0.5)), mytable.Tremble_Dur(mytable.Ingress_Mag<0.5), 100, 'filled', 'markerfacealpha',0.5 ); set(gca,'Colormap',[0,0,0;1,0,0]); xlim([-10,5]); xlabel('TrembleMag');ylabel('TrembleDur'); title('No ingress'); ylim([0,1.1])

%figure('color','w');
%mytable = summarize( mymice );
%subplot(1,2,1); scatter( mytable.Tremble_Mag(mytable.Ingress_Mag>0.5)), mytable.Tremble_Dur(mytable.Ingress_Mag>0.5), 100, 'filled', 'markerfacealpha',0.5 ); set(gca,'Colormap',[0,0,0;1,0,0]); xlim([-10,5]); xlabel('TrembleMag');ylabel('TrembleDur'); title('Ingress'); ylim([0,1.1])
%subplot(1,2,2); scatter( mytable.Tremble_Mag(mytable.Ingress_Mag<0.5)), mytable.Tremble_Dur(mytable.Ingress_Mag<0.5), 100, 'filled', 'markerfacealpha',0.5 ); set(gca,'Colormap',[0,0,0;1,0,0]); xlim([-10,5]); xlabel('TrembleMag');ylabel('TrembleDur'); title('No ingress'); ylim([0,1.1])


%%

figure('color','w');

mytable = summarize( mymice )
mytable = mytable( or(mytable.TrialType=='CS+1',mytable.TrialType=='CS+2'), : )
subplot(2,2,1); scatter( log(mytable.Tremble_Mag(mytable.Ingress_Mag>0.5)), mytable.Tremble_Dur(mytable.Ingress_Mag>0.5), 100, mytable.Trials(mytable.Ingress_Mag>0.5), 'filled', 'markerfacealpha',0.5, 'markeredgecolor', 'k' ); set(gca,'Colormap',gray(10)); xlim([-10,5]); xlabel('TrembleMag');ylabel('TrembleDur'); title('Ingress (CS+)'); ylim([0,1.1])
subplot(2,2,2); scatter( log(mytable.Tremble_Mag(mytable.Ingress_Mag<0.5)), mytable.Tremble_Dur(mytable.Ingress_Mag<0.5), 100, mytable.Trials(mytable.Ingress_Mag<0.5), 'filled', 'markerfacealpha',0.5, 'markeredgecolor', 'k' ); set(gca,'Colormap',gray(10)); xlim([-10,5]); xlabel('TrembleMag');ylabel('TrembleDur'); title('No Ingress (CS+)'); ylim([0,1.1])

mytable = summarize( mymice )
mytable = mytable( or(mytable.TrialType=='CS-1',mytable.TrialType=='CS-2'), : )
subplot(2,2,3); scatter( log(mytable.Tremble_Mag(mytable.Ingress_Mag>0.5)), mytable.Tremble_Dur(mytable.Ingress_Mag>0.5), 100, mytable.Trials(mytable.Ingress_Mag>0.5), 'filled', 'markerfacealpha',0.5, 'markeredgecolor', 'k' ); set(gca,'Colormap',gray(10)); xlim([-10,5]); xlabel('TrembleMag');ylabel('TrembleDur'); title('Ingress (CS-)'); ylim([0,1.1])
subplot(2,2,4); scatter( log(mytable.Tremble_Mag(mytable.Ingress_Mag<0.5)), mytable.Tremble_Dur(mytable.Ingress_Mag<0.5), 100, mytable.Trials(mytable.Ingress_Mag<0.5), 'filled', 'markerfacealpha',0.5, 'markeredgecolor', 'k' ); set(gca,'Colormap',gray(10)); xlim([-10,5]); xlabel('TrembleMag');ylabel('TrembleDur'); title('No Ingress (CS-)'); ylim([0,1.1])
%suptitle('CS-')

%%

figure('color','w');

mytable = summarize( mymice )
mytable = mytable( and( mytable.Ingress_onset>0.2, or(mytable.TrialType=='CS+1',mytable.TrialType=='CS+2')), : )
subplot(2,2,1); scatter( log(mytable.Tremble_Mag(mytable.Ingress_Mag>0.7)), mytable.Tremble_Dur(mytable.Ingress_Mag>0.7), 100*mytable.Ingress_onset(mytable.Ingress_Mag>0.7), mytable.Mice(mytable.Ingress_Mag>0.7), 'filled', 'markerfacealpha',0.7, 'markeredgecolor', 'k' ); set(gca,'Colormap',jet(20)); xlim([-10,5]); xlabel('TrembleMag');ylabel('TrembleDur'); title('Ingress (CS+)'); ylim([0,1.1])
subplot(2,2,2); scatter( log(mytable.Tremble_Mag(mytable.Ingress_Mag<0.7)), mytable.Tremble_Dur(mytable.Ingress_Mag<0.7), 100*mytable.Ingress_onset(mytable.Ingress_Mag<0.7), mytable.Mice(mytable.Ingress_Mag<0.7), 'filled', 'markerfacealpha',0.7, 'markeredgecolor', 'k' ); set(gca,'Colormap',jet(20)); xlim([-10,5]); xlabel('TrembleMag');ylabel('TrembleDur'); title('No Ingress (CS+)'); ylim([0,1.1])

mytable = summarize( mymice )
mytable = mytable( and( mytable.Ingress_onset>0.2, or(mytable.TrialType=='CS-1',mytable.TrialType=='CS-2')), : )
subplot(2,2,3); scatter( log(mytable.Tremble_Mag(mytable.Ingress_Mag>0.7)), mytable.Tremble_Dur(mytable.Ingress_Mag>0.7), 100*mytable.Ingress_onset(mytable.Ingress_Mag>0.7), mytable.Mice(mytable.Ingress_Mag>0.7), 'filled', 'markerfacealpha',0.7, 'markeredgecolor', 'k' ); set(gca,'Colormap',jet(20)); xlim([-10,5]); xlabel('TrembleMag');ylabel('TrembleDur'); title('Ingress (CS-)'); ylim([0,1.1])
subplot(2,2,4); scatter( log(mytable.Tremble_Mag(mytable.Ingress_Mag<0.7)), mytable.Tremble_Dur(mytable.Ingress_Mag<0.7), 100*mytable.Ingress_onset(mytable.Ingress_Mag<0.7), mytable.Mice(mytable.Ingress_Mag<0.7), 'filled', 'markerfacealpha',0.7, 'markeredgecolor', 'k' ); set(gca,'Colormap',jet(20)); xlim([-10,5]); xlabel('TrembleMag');ylabel('TrembleDur'); title('No Ingress (CS-)'); ylim([0,1.1])
%suptitle('CS-')

%% Tremble mag versus ingress mag by mouse

mytable = summarize( mymice )
condition = or(mytable.TrialType=='CS+1',mytable.TrialType=='CS+2');
summary = varfun( @nanmean, mytable(condition,{'Mice','Tremble_Mag','Ingress_Mag','Ingress_onset'}), 'GroupingVar', 'Mice' );

figure('color','w'); scatter( summary.nanmean_Tremble_Mag, summary.nanmean_Ingress_Mag, 100, 'filled', 'markeredgecolor','k' );
rowfun( @(mouse,~,x,y,~) text( x, y, sprintf('%s',mouse) ), summary )
set(gca,'TickDir','out'); xlabel('Tremble mag'); ylabel('Ingress mag')

%% Tremble mag versus ingress onset by mouse

mytable = summarize( mymice )
condition = or(mytable.TrialType=='CS+1',mytable.TrialType=='CS+2');
summary = varfun( @nanmedian, mytable(condition,{'Mice','Tremble_Mag','Ingress_onset','Tremble_Dur'}), 'GroupingVar', 'Mice' );

figure('color','w'); scatter( summary.nanmedian_Tremble_Mag, summary.nanmedian_Ingress_onset, 100, summary.nanmedian_Tremble_Dur, 'filled', 'markeredgecolor','k' );
rowfun( @(mouse,~,x,y,~) text( x, y, sprintf('%s',mouse) ), summary )
set(gca,'TickDir','out'); xlabel('Tremble Mag'); ylabel('Ingress Onset')


%% Ingress onset by mouse

mytable = summarize( mymice )
condition = or(mytable.TrialType=='CS+1',mytable.TrialType=='CS+2');
mytable = mytable( condition, : );
figure; boxplot( mytable.Ingress_onset, mytable.Mice )

figure('color','w'); scatter( summary.nanmean_Ingress_onset, summary.nanmean_Tremble_Mag, 100, 'filled', 'markeredgecolor','k' );
set(gca,'TickDir','out'); xlabel('Ingress onset'); ylabel('Tremble Mag')

%%

figure('color','w'); scatter( summary.nanmean_Ingress_Mag, summary.nanmean_Tremble_Mag, 100, 'filled', 'markeredgecolor','k' );
set(gca,'TickDir','out'); xlabel('Ingress Mag'); ylabel('Tremble Mag')

%%


mytable = summarize( mymice )
mytable = mytable( or(mytable.TrialType=='CS+1',mytable.TrialType=='CS+2'), : )
figure('color','w'); scatter( log( mytable.Tremble_Mag ), mytable.Ingress_Mag )
set(gca,'TickDir','out'); ylabel('Ingress Mag'); xlabel('Tremble Mag')

mytable = summarize( mymice )
mytable = mytable( or(mytable.TrialType=='CS-1',mytable.TrialType=='CS-2'), : )
figure('color','w'); scatter( log( mytable.Tremble_Mag ), mytable.Ingress_Mag )
set(gca,'TickDir','out'); ylabel('Ingress Mag'); xlabel('Tremble Mag')


%% DREADD Animals


mytable = mytable_dreadd
mytable.TrialType = categorical( mytable.TrialType )
mytable.Mice = categorical( mytable.Mice )
condition = or(mytable.TrialType=='CS+1',mytable.TrialType=='CS+2');
summary = varfun( @nanmean, mytable(condition,{'Mice','Tremble_Mag','Ingress_Mag','Ingress_onset'}), 'GroupingVar', 'Mice' );

figure('color','w'); scatter( summary.nanmean_Tremble_Mag, summary.nanmean_Ingress_Mag, 100, 'filled', 'markeredgecolor','k' );
rowfun( @(mouse,~,x,y,~) text( x, y, sprintf('%s',mouse) ), summary )
set(gca,'TickDir','out'); xlabel('Tremble mag'); ylabel('Ingress mag')