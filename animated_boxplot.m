%% CS+1 boxplot animation

figure('color','w')
trialType = 'CS+1';

mytable = summarize( mymice );
mytable = mytable( mytable.TrialType==trialType, : )
boxplot( mytable.Ingress_Mag, mytable.Trials ); hold on;
scatter( mytable.Trials, mytable.Ingress_Mag, 50, mytable.Mice, 'filled', 'markerfacealpha', 0.5, 'linewidth', 1, 'markeredgecolor', 'none' )

boxes = findobj('tag','Box'); arrayfun( @(x) set(x,'Color','k'), boxes )
boxes = findobj('tag','Median'); arrayfun( @(x) set(x,'Color','k'), boxes )

box off; set(gca,'TickDir','out')
i=1;

for mice_ = unique( mytable.Mice )'
    
    mytable_ = mytable( and( mytable.Mice==mice_,mytable.TrialType==trialType ), : );
    t1 = scatter( mytable_.Trials, mytable_.Ingress_Mag, 100, mytable_.Mice, 'filled', 'markerfacecolor', 'r', 'markerfacealpha', 1, 'linewidth', 1, 'markeredgecolor', 'k' );
    t2 = text( 10.5, mytable_.Ingress_Mag(end), mice_ );
    ylim([0,1]);
    
    filename = 'animated_extinction_csplus1_bp.gif';
    
    drawnow 
    % Capture the plot as an image 
    frame = getframe(gcf); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    xlabel('Trial'); ylabel('Ingress magnitude');
    % Write to the GIF File 
    if i == 1 
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      fprintf('%i\n',i)
    end
    xlabel('Trial'); ylabel('Ingress magnitude');
    
    i=i+1;
    pause(1);
    delete(t1); delete(t2);
end

% CS-1 boxplot animation

figure('color','w')
trialType = 'CS-1';

mytable = summarize( mymice );
mytable = mytable( mytable.TrialType==trialType, : )
boxplot( mytable.Ingress_Mag, mytable.Trials ); hold on;
scatter( mytable.Trials, mytable.Ingress_Mag, 50, mytable.Mice, 'filled', 'markerfacealpha', 0.5, 'linewidth', 1, 'markeredgecolor', 'none' )

boxes = findobj('tag','Box'); arrayfun( @(x) set(x,'Color','k'), boxes )
boxes = findobj('tag','Median'); arrayfun( @(x) set(x,'Color','k'), boxes )

box off; set(gca,'TickDir','out')
i=1;

for mice_ = unique( mytable.Mice )'
    
    mytable_ = mytable( and( mytable.Mice==mice_,mytable.TrialType==trialType ), : );
    t1 = scatter( mytable_.Trials, mytable_.Ingress_Mag, 100, mytable_.Mice, 'filled', 'markerfacecolor', 'r', 'markerfacealpha', 1, 'linewidth', 1, 'markeredgecolor', 'k' )
    t2 = text( 10.5, mytable_.Ingress_Mag(end), mice_ )
    ylim([0,1]);
    
    filename = 'animated_extinction_csminus1_bp.gif';
    
    drawnow 
    % Capture the plot as an image 
    frame = getframe(gcf); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    xlabel('Trial'); ylabel('Ingress magnitude');
    % Write to the GIF File 
    if i == 1 
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      fprintf('%i\n',i)
    end
    
    i=i+1;
    pause(1);
    delete(t1); delete(t2);
end

%% CS+1 boxplot animation

figure('color','w')
trialType = 'CS+1';

mytable = summarize( mymice );
mytable = mytable( mytable.TrialType==trialType, : )
boxplot( mytable.Tremble_Dur, mytable.Trials ); hold on;
scatter( mytable.Trials, mytable.Tremble_Dur, 50, mytable.Mice, 'filled', 'markerfacealpha', 0.5, 'linewidth', 1, 'markeredgecolor', 'none' )

boxes = findobj('tag','Box'); arrayfun( @(x) set(x,'Color','k'), boxes )
boxes = findobj('tag','Median'); arrayfun( @(x) set(x,'Color','k'), boxes )

box off; set(gca,'TickDir','out')
i=1;

for mice_ = unique( mytable.Mice )'
    
    mytable_ = mytable( and( mytable.Mice==mice_,mytable.TrialType==trialType ), : );
    t1 = scatter( mytable_.Trials, mytable_.Tremble_Dur, 100, mytable_.Mice, 'filled', 'markerfacecolor', 'r', 'markerfacealpha', 1, 'linewidth', 1, 'markeredgecolor', 'k' );
    t2 = text( 10.5, mytable_.Tremble_Dur(end), mice_ );
    ylim([0,1]);
    
    filename = 'animated_tremble_extinction_csplus1_bp.gif';
    
    drawnow 
    % Capture the plot as an image 
    frame = getframe(gcf); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    xlabel('Trial'); ylabel('Ingress magnitude');
    % Write to the GIF File 
    if i == 1 
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      fprintf('%i\n',i)
    end
    xlabel('Trial'); ylabel('Tremble magnitude');
    
    i=i+1;
    pause(1);
    delete(t1); delete(t2);
end

% CS-1 boxplot animation

figure('color','w')
trialType = 'CS-1';

mytable = summarize( mymice );
mytable = mytable( mytable.TrialType==trialType, : )
boxplot( mytable.Tremble_Dur, mytable.Trials ); hold on;
scatter( mytable.Trials, mytable.Tremble_Dur, 50, mytable.Mice, 'filled', 'markerfacealpha', 0.5, 'linewidth', 1, 'markeredgecolor', 'none' )

boxes = findobj('tag','Box'); arrayfun( @(x) set(x,'Color','k'), boxes )
boxes = findobj('tag','Median'); arrayfun( @(x) set(x,'Color','k'), boxes )

box off; set(gca,'TickDir','out')
i=1;

for mice_ = unique( mytable.Mice )'
    
    mytable_ = mytable( and( mytable.Mice==mice_,mytable.TrialType==trialType ), : );
    t1 = scatter( mytable_.Trials, mytable_.Tremble_Dur, 100, mytable_.Mice, 'filled', 'markerfacecolor', 'r', 'markerfacealpha', 1, 'linewidth', 1, 'markeredgecolor', 'k' )
    t2 = text( 10.5, mytable_.Tremble_Dur(end), mice_ )
    ylim([0,1]);
    
    filename = 'animated_tremble_extinction_csminus1_bp.gif';
    
    drawnow 
    % Capture the plot as an image 
    frame = getframe(gcf); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    xlabel('Trial'); ylabel('Tremble magnitude');
    % Write to the GIF File 
    if i == 1 
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      fprintf('%i\n',i)
    end
    
    i=i+1;
    pause(1);
    delete(t1); delete(t2);
end

%%

% All data
mytable = summarize( mymice );
%mytable = mytable( or( mytable.TrialType=='CS+1',mytable.TrialType=='CS+2' ), : )

output = varfun( @nanmean, mytable(:,[1,5,7]), 'GroupingVariables', {'Mice'} )

figure('color','w');
scatter( output.nanmean_Ingress_Mag, output.nanmean_Tremble_Dur, 100,'filled', 'markeredgecolor', 'k' )
set(gca,'TickDir','out','xlim',[0,1],'ylim',[0,1],'xlabel',text(0,0,'Ingress Magnitude'),'ylabel',text(0,0,'Tremble Duration') )
line([0,1],[0,1],'color','k')

%%

% All data
mytable = summarize( mymice );
%mytable = mytable( or( mytable.TrialType=='CS+1',mytable.TrialType=='CS+2' ), : )


figure('color','w');
scatter( mytable.Ingress_Mag, mytable.Tremble_Dur, 100, mytable.Mice, 'filled', 'markeredgecolor', 'k', 'markerfacealpha', 0.5 )
set(gca,'TickDir','out','xlim',[0,1],'ylim',[0,1],'xlabel',text(0,0,'Ingress Magnitude'),'ylabel',text(0,0,'Tremble Duration') )
line([0,1],[0,1],'color','k')
