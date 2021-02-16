%% Logistic regression tutorial
rng(1)
low_tremble_scores = min(.3,normrnd(.3,.1,20,1));
high_tremble_scores = min(1,normrnd(.5,.25,20,1));
low_breathing_scores = min(.3,normrnd(.3,.1,20,1));
high_breathing_scores = min(1,normrnd(.5,.25,20,1));

tremble_scores = [low_tremble_scores;high_tremble_scores];
breathing_scores = [low_breathing_scores;high_breathing_scores];
ingress = [zeros(19,1);1;ones(19,1);0];

mouse_tbl = table( tremble_scores,breathing_scores,ingress )
mdl = fitglm( mouse_tbl, modelspec, 'Distribution', 'binomial' )


%%
figure()
ax = axes('nextplot','add','tickdir','out')
all_fp = [];
all_tp = [];

for threshold = [0:.01:1]
   tp = intersect( find( OR./(1+OR) > threshold ) , find( mouse_tbl.ingress == 1 ) );
   fp = intersect( find( OR./(1+OR) > threshold ) , find( mouse_tbl.ingress == 0 ) );
   fprintf('TP: %i FP: %i\n',numel(tp),numel(fp))
   all_fp = [ all_fp; numel(fp) ];
   all_tp = [ all_tp; numel(tp) ];
   plot( numel(fp), numel(tp), 'o' )
end

[all_fp,sort_all_fp] = sort( all_fp );
figure; plot( all_fp, all_tp( sort_all_fp ) );

AUROC = cumtrapz( all_fp,  all_tp( sort_all_fp ) )./400;
AUROC = AUROC(end)
%% Logistic regression only for CS- trials
mytable_ = summarize(mymice);
modelspec = 'Ingress_Mag ~ Tremble_Mag + Tremble_Dur'
mytable_.Ingress_Mag = mytable_( :,'Ingress_Mag').Ingress_Mag>.5

tmptable = mytable_(mytable_.TrialType_Generic=='CS+',:);
mdl = fitglm( tmptable, modelspec, 'Distribution', 'binomial' )

%%

[y1,x1] = ecdf( mytable_(mytable_.TrialType_Generic=='CS+',:).Tremble_Dur );
[y2,x2] = ecdf( mytable_(mytable_.TrialType_Generic=='CS-',:).Tremble_Dur );

figure; plot(x1,y1); hold on; plot(x2,y2);

%% FULL model

mytable_ = summarize(mymice);
mytable__ = mytable_( :, : );
mytable__.Ingress_Mag = mytable__(:,'Ingress_Mag').Ingress_Mag>.5
modelspec = 'Ingress_Mag ~ Tremble_Mag + Tremble_Dur + Ingress_onset'
mdl = fitglm( mytable__, modelspec, 'Distribution', 'binomial' )
mylogit = @(x1,x2,x3) exp( -0.50466-1.4161*x1+2.057*x2-0.76443*x3 )

output = rowfun( @(x1,x2,x3) mylogit(x1,x2,x3), mytable__(:,{'Tremble_Mag','Tremble_Dur','Ingress_onset'}) );
OR = output.Var1;

all_fp = [];
all_tp = [];

mytable__.OR = OR;

for threshold = [0:.001:1]
   tp = intersect( find( mytable__.OR./(1+mytable__.OR) > threshold ) , find( mytable__.Ingress_Mag == 1 ) );
   fp = intersect( find( mytable__.OR./(1+mytable__.OR) > threshold ) , find( mytable__.Ingress_Mag == 0 ) );
   fn = intersect( find( mytable__.OR./(1+mytable__.OR) < threshold ) , find( mytable__.Ingress_Mag == 1 ) );
   tn = intersect( find( mytable__.OR./(1+mytable__.OR) < threshold ) , find( mytable__.Ingress_Mag == 0 ) );
   fprintf('TP: %i FP: %i\n',numel(tp),numel(fp))
   all_fp = [ all_fp; numel(fp)./( numel(fp)+numel(tn) ) ];
   all_tp = [ all_tp; numel(tp)./( numel(tp)+numel(fn) ) ];
end

[all_fp,sort_all_fp] = sort( all_fp );
figure; plot( all_fp, all_tp( sort_all_fp ) );

AUROC = cumtrapz( all_fp,  all_tp( sort_all_fp ) )
AUROC = AUROC(end)


%% CS- model

mytable_ = summarize(mymice);
mytable__ = mytable_( mytable_.TrialType_Generic=='CS-', : );
modelspec = 'Ingress_Mag ~ Tremble_Mag + Tremble_Dur'
mdl = fitglm( mytable__, modelspec, 'Distribution', 'binomial' )
mylogit = @(x1,x2) exp( -2.4685+1.4513*x1+1.1247*x2 )
mytable__.Ingress_Mag = mytable__(:,'Ingress_Mag').Ingress_Mag>.5

[y1,x1] = ecdf( mytable__( mytable__.Ingress_Mag==1, 'Tremble_Dur' ).Tremble_Dur );
[y2,x2] = ecdf( mytable__( mytable__.Ingress_Mag==0, 'Tremble_Dur' ).Tremble_Dur );
figure; plot( x1,y1 ); hold on; plot( x2, y2 ); 
legend({'Ingress','No Ingress'}); xlabel('Tremble Dur'); ylabel('Fraction of trials')

output = rowfun( @(x1,x2) mylogit(x1,x2), mytable__(:,{'Tremble_Mag','Tremble_Dur'}) );
OR = output.Var1;

all_fp = [];
all_tp = [];

mytable__.OR = OR;

for threshold = [0:.001:1]
   tp = intersect( find( mytable__.OR./(1+mytable__.OR) > threshold ) , find( mytable__.Ingress_Mag == 1 ) );
   fp = intersect( find( mytable__.OR./(1+mytable__.OR) > threshold ) , find( mytable__.Ingress_Mag == 0 ) );
   fn = intersect( find( mytable__.OR./(1+mytable__.OR) < threshold ) , find( mytable__.Ingress_Mag == 1 ) );
   tn = intersect( find( mytable__.OR./(1+mytable__.OR) < threshold ) , find( mytable__.Ingress_Mag == 0 ) );
   fprintf('TP: %i FP: %i\n',numel(tp),numel(fp))
   all_fp = [ all_fp; numel(fp)./( numel(fp)+numel(tn) ) ];
   all_tp = [ all_tp; numel(tp)./( numel(tp)+numel(fn) ) ];
end

[all_fp,sort_all_fp] = sort( all_fp );
figure; plot( all_fp, all_tp( sort_all_fp ) );

AUROC = cumtrapz( all_fp,  all_tp( sort_all_fp ) )
AUROC = AUROC(end)


%% CS+ model

mytable_ = summarize(mymice);
mytable__ = mytable_( mytable_.TrialType_Generic=='CS+', : );
mytable__.Ingress_Mag = mytable__(:,'Ingress_Mag').Ingress_Mag>.5
modelspec = 'Ingress_Mag ~ Tremble_Mag + Tremble_Dur'
mdl = fitglm( mytable__, modelspec, 'Distribution', 'binomial' )
mylogit = @(x1,x2) exp( -0.91497-2.5507*x1+2.6136*x2 )

%
[y1,x1] = ecdf( mytable__( mytable__.Ingress_Mag==1, 'Tremble_Dur' ).Tremble_Dur );
[y2,x2] = ecdf( mytable__( mytable__.Ingress_Mag==0, 'Tremble_Dur' ).Tremble_Dur );
figure; plot( x1,y1 ); hold on; plot( x2, y2 ); 
legend({'Ingress','No Ingress'}); xlabel('Tremble Dur'); ylabel('Fraction of trials')

output = rowfun( @(x1,x2) mylogit(x1,x2), mytable__(:,{'Tremble_Mag','Tremble_Dur'}) );
OR = output.Var1;

all_fp = [];
all_tp = [];

mytable__.OR = OR;

for threshold = [0:.001:1]
   tp = intersect( find( mytable__.OR./(1+mytable__.OR) > threshold ) , find( mytable__.Ingress_Mag == 1 ) );
   fp = intersect( find( mytable__.OR./(1+mytable__.OR) > threshold ) , find( mytable__.Ingress_Mag == 0 ) );
   fn = intersect( find( mytable__.OR./(1+mytable__.OR) < threshold ) , find( mytable__.Ingress_Mag == 1 ) );
   tn = intersect( find( mytable__.OR./(1+mytable__.OR) < threshold ) , find( mytable__.Ingress_Mag == 0 ) );
   fprintf('TP: %i FP: %i\n',numel(tp),numel(fp))
   all_fp = [ all_fp; numel(fp)./( numel(fp)+numel(tn) ) ];
   all_tp = [ all_tp; numel(tp)./( numel(tp)+numel(fn) ) ];
end

[all_fp,sort_all_fp] = sort( all_fp );
figure; plot( all_fp, all_tp( sort_all_fp ) );

AUROC = cumtrapz( all_fp,  all_tp( sort_all_fp ) )
AUROC = AUROC(end)

%% Trying it on the real data
mytable = summarize(mymice);
mytable_ = mytable(:,{'mouseIds','TrialType_Generic','Ingress_Mag','Tremble_Mag','Tremble_Dur'});
mytable_.Ingress_Mag = mytable_(:,'Ingress_Mag').Ingress_Mag>.5
modelspec = 'Ingress_Mag ~ Tremble_Mag + Tremble_Dur'

mdl = fitglm( mytable_, modelspec, 'Distribution', 'binomial' )

mylogit = @(x1,x2) exp(-1.7689 -6.7003*x1+3.6053*x2 )

output = rowfun( @(x1,x2) mylogit(x1,x2), mytable_(:,{'Tremble_Mag','Tremble_Dur'}) );
OR = output.Var1;


figure()
ax = axes('nextplot','add','tickdir','out')

%% CS type comparison

all_fp = [];
all_tp = [];

mytable_.OR = OR;
tmptable = mytable_(mytable_.TrialType_Generic=='CS-',:);

for threshold = [0:.01:1]
   tp = intersect( find( tmptable.OR./(1+tmptable.OR) > threshold ) , find( tmptable.Ingress_Mag == 1 ) );
   fp = intersect( find( tmptable.OR./(1+tmptable.OR) > threshold ) , find( tmptable.Ingress_Mag == 0 ) );
   fn = intersect( find( tmptable.OR./(1+tmptable.OR) < threshold ) , find( tmptable.Ingress_Mag == 1 ) );
   tn = intersect( find( tmptable.OR./(1+tmptable.OR) < threshold ) , find( tmptable.Ingress_Mag == 0 ) );
   fprintf('TP: %i FP: %i\n',numel(tp),numel(fp))
   all_fp = [ all_fp; numel(fp)./( numel(fp)+numel(tn) ) ];
   all_tp = [ all_tp; numel(tp)./( numel(tp)+numel(fn) ) ];
end

[all_fp,sort_all_fp] = sort( all_fp );
figure; plot( all_fp, all_tp( sort_all_fp ) );

AUROC = cumtrapz( all_fp,  all_tp( sort_all_fp ) )
AUROC = AUROC(end)



%%

all_fp = [];
all_tp = [];

mouseKey = 'PO371';

function calc_AUROC( mytable_, mouseKey )
mytable_.OR = OR;
tmptable = mytable_(mytable_.mouseIds==mousekey,:);

for threshold = [0:.01:1]
   tp = intersect( find( tmptable.OR./(1+tmptable.OR) > threshold ) , find( tmptable.Ingress_Mag == 1 ) );
   fp = intersect( find( tmptable.OR./(1+tmptable.OR) > threshold ) , find( tmptable.Ingress_Mag == 0 ) );
   fn = intersect( find( tmptable.OR./(1+tmptable.OR) < threshold ) , find( tmptable.Ingress_Mag == 1 ) );
   tn = intersect( find( tmptable.OR./(1+tmptable.OR) < threshold ) , find( tmptable.Ingress_Mag == 0 ) );
   fprintf('TP: %i FP: %i\n',numel(tp),numel(fp))
   all_fp = [ all_fp; numel(fp)./( numel(fp)+numel(tn) ) ];
   all_tp = [ all_tp; numel(tp)./( numel(tp)+numel(fn) ) ];
end

[all_fp,sort_all_fp] = sort( all_fp );
figure; plot( all_fp, all_tp( sort_all_fp ) );

AUROC = cumtrapz( all_fp,  all_tp( sort_all_fp ) )
AUROC = AUROC(end)
end
