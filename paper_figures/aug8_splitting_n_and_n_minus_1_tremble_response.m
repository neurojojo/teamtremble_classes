mouse = 'PO352';

keystr = mymice.mouseObjs.(mouse).keyString;

figure('color','w','windowstate','fullscreen');
minus_plus = regexp( keystr, '[pm][12]', 'match' );
find( arrayfun( @(x) strcmp(x,'m1'), minus_plus ) == 1 );
find( arrayfun( @(x) strcmp(x,'m2'), minus_plus ) == 1 );
find( arrayfun( @(x) strcmp(x,'p1'), minus_plus ) == 1 );
find( arrayfun( @(x) strcmp(x,'p2'), minus_plus ) == 1 );

newtable( find( arrayfun( @(x) strcmp(x,'m2'), minus_plus ) == 1 ), : ) = mytable( and(mytable.Mice==mouse,mytable.TrialType=='CS-2'), : );
newtable( find( arrayfun( @(x) strcmp(x,'m1'), minus_plus ) == 1 ), : ) = mytable( and(mytable.Mice==mouse,mytable.TrialType=='CS-1'), : );
newtable( find( arrayfun( @(x) strcmp(x,'p1'), minus_plus ) == 1 ), : ) = mytable( and(mytable.Mice==mouse,mytable.TrialType=='CS+1'), : );
newtable( find( arrayfun( @(x) strcmp(x,'p2'), minus_plus ) == 1 ), : ) = mytable( and(mytable.Mice==mouse,mytable.TrialType=='CS+2'), : );

scatter( newtable.Tremble_Dur(1:end-1), newtable.Tremble_Dur(2:end), 150, newtable.TrialType(2:end), 'filled' ); hold on; scatter( newtable.Tremble_Dur(1:end-1), newtable.Tremble_Dur(2:end), 50, newtable.TrialType(1:end-1), 'filled' ); xlim([0,1]); ylim([0,1])
set(gca,'ColorMap',[0,0.2,1;0,0,1;.8,.8,.8;.7,.7,.7],'TickDir','out','FontSize',16,'xlim',[0,1.1],'ylim',[0,1.1]);
plot([0,1],[0,1],'k--','linewidth',4);
title(mouse,'FontWeight','normal')

arrayfun( @(x) text( newtable(x,:).Tremble_Dur+.005, newtable(x+1,:).Tremble_Dur+.005, sprintf('%i',x+1) ), [1:39] )


figure; bar( [1:40], newtable.Tremble_Dur ); set(gca,'XTick',[1:40],'XTickLabel', newtable.TrialType, 'XTickLabelRotation',90)


%% Previous 2 versus current 

results = structfun( @(x) x.searchMouse('p[12]p[12]m[12]'), mymice.mouseObjs, 'UniformOutput', false )

figure;
ax1 = axes('nextplot','add');

for mouse = fields(results)'
   
    array_ = table2array(results.(mouse{1}));
    if isempty(array_); continue; end;
    data_ = [ mymice.mouseObjs.(mouse{1}).trembleScore_csplus.Dur( array_(:,[1:2]) ), mymice.mouseObjs.(mouse{1}).trembleScore_csminus.Dur( array_(:,3) )' ]
    plot( ax1, nanmean(data_,1) );
    
end

%% Plus, minus, minus, plus


results = structfun( @(x) x.searchMouse('p[12]m[12]m[12]p[12]'), mymice.mouseObjs, 'UniformOutput', false )

figure('color','w');
ax1 = axes('nextplot','add','ylim',[0,1],'tickdir','out','xlim',[0,5],'xtick',[1:4],'xticklabel',{'CS+','CS-','CS-','CS+'});
data_all = [];

for mouse = fields(results)'
   
    array_ = table2array(results.(mouse{1}));
    if isempty(array_); continue; end;
    data_ = [ mymice.mouseObjs.(mouse{1}).trembleScore_csplus.Dur( array_(:,[1]) )', mymice.mouseObjs.(mouse{1}).trembleScore_csminus.Dur( array_(:,[2:3]) ), mymice.mouseObjs.(mouse{1}).trembleScore_csplus.Dur( array_(:,[4]) )' ];
    plot( ax1, [1:4], nanmean(data_,1), '.-', 'markersize', 24 );
    data_all = [ nanmean(data_,1); data_all ];
    
end
ylabel('Tremble Duration')
saveas(gcf,'Tremble_Duration_csp_csm_csm_csp.svg')

%% Plus, minus, minus, plus


results = structfun( @(x) x.searchMouse('p[12]m[12]m[12]p[12]'), mymice.mouseObjs, 'UniformOutput', false )

figure('color','w');
ax1 = axes('nextplot','add','ylim',[0,1],'tickdir','out','xlim',[0,5],'xtick',[1:4],'xticklabel',{'CS+','CS-','CS-','CS+'});
data_all = [];

for mouse = fields(results)'
   
    array_ = table2array(results.(mouse{1}));
    if isempty(array_); continue; end;
    data_ = [ mymice.mouseObjs.(mouse{1}).trembleScore_csplus.Dur( array_(:,[1]) )', mymice.mouseObjs.(mouse{1}).trembleScore_csminus.Dur( array_(:,[2:3]) ), mymice.mouseObjs.(mouse{1}).trembleScore_csplus.Dur( array_(:,[4]) )' ];
    plot( ax1, [1:4], nanmean(data_,1), '.-', 'markersize', 24 );
    data_all = [ nanmean(data_,1); data_all ];
    
end
ylabel('Tremble Duration')
saveas(gcf,'Tremble_Duration_csp_csm_csm_csp.svg')

%% Looking at the averages of ingress magnitude and tremble duration

output = varfun( @nanmean, mytable( or(mytable.TrialType=='CS+1',mytable.TrialType=='CS+2'),[1,5,6,7]), 'GroupingVar', 'Mice' )
figure; scatter( output.nanmean_Ingress_Mag, output.nanmean_Tremble_Dur, min(5,output.nanmean_Tremble_Mag)*50, 'filled'  );

rowfun( @(mouse,grp,mag,tmag,dur) text( mag, dur, mouse ), output )


%%

output = varfun( @nanmean, mytable( :,[1,5,6,7,9]), 'GroupingVar', {'Mice','TrialType_Generic'} )

figure;
scatter( output.nanmean_Ingress_Mag, output.nanmean_Tremble_Dur, max( 1, min(5,output.nanmean_Tremble_Mag) )*50, output.TrialType_Generic, 'filled'  );
xlabel('IngressMag'); ylabel('TrembleDur'); title('Size = Tremble Mag')

rowfun( @(mouse,type,grp,mag,tmag,dur) text( mag+0.01, dur+0.01, mouse ), output )
