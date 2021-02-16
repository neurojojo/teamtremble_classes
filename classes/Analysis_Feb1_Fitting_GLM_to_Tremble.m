%%

Ndata = size( mytable, 1 );
train_set = randperm(Ndata,540);
held_out = setdiff( [1:Ndata], train_set );
tmp = mytable( mytable.Tremble_Dur>.2, : );
mdl = fitglm( tmp, modelspec, 'Distribution', 'normal' )

%%
tmp = mytable( mytable.Tremble_Dur>.2, : );
Nheld_out = size(tmp,1);
w = [ 0.25759, -0.69254, 0.44782, 0.016634]';

predicted_ingress_mag = [ones(Nheld_out,1),tmp.Tremble_Mag,tmp.Tremble_Dur,tmp.Ingress_onset]*w;

figure; scatter( tmp.Ingress_Mag, predicted_ingress_mag, 100, tmp.mouseIds, 'filled' );
xlabel('Actual'); ylabel('Predicted')
xlim([0,1]);ylim([0,1]);line([0,1],[0,1])

tmp.predicted_ingress_mag = predicted_ingress_mag
tmp.error = abs( tmp.Ingress_Mag - predicted_ingress_mag )

%tmp = tmp( :, {'mouseIds','TrialType','TrialType_Generic','TrialNumber','Ingress_Mag','predicted_ingress_mag','Tremble_Mag','Tremble_Dur','Ingress_onset','Ingress','error'})

%%
rng(1)
Ndata = size( mytable, 1 );
tmp = mytable;
mdl = fitglm( tmp, modelspec, 'Distribution', 'normal' )

%%
tmp = mytable;
tmp = tmp( isnan(tmp.Ingress_onset)==0, : );
Nheld_out = size(tmp,1);
w = [ 0.42116, -0.45288, 0.26658, -0.14854]';

predicted_ingress_mag = [ones(Nheld_out,1),tmp.Tremble_Mag,tmp.Tremble_Dur,tmp.Ingress_onset]*w;

figure; scatter( tmp.Ingress_Mag(tmp.TrialNumber==5), predicted_ingress_mag(tmp.TrialNumber==5), 100, tmp.mouseIds(tmp.TrialNumber==5), 'filled' );
xlabel('Actual'); ylabel('Predicted')
xlim([0,1]);ylim([0,1]);line([0,1],[0,1])

tmp.predicted_ingress_mag = predicted_ingress_mag
tmp.error = abs( tmp.Ingress_Mag - predicted_ingress_mag )

tmp = tmp( :, {'mouseIds','TrialType','TrialType_Generic','TrialNumber','Ingress_Mag','predicted_ingress_mag','Tremble_Mag','Tremble_Dur','Ingress_onset','Ingress','error'})
%%
figure;boxplot( tmp.error, tmp.TrialNumber ); figure; scatter( tmp.TrialNumber, tmp.error, 100, tmp.TrialType_Generic, 'filled', 'jitter', 'on' )
figure; scatter( tmp.TrialNumber, tmp.error, 100, tmp.mouseIds, 'filled', 'jitter', 'on' )
%% Shuffled

rng(1)
Ndata = size( mytable, 1 );
tmp = mytable;
tmp.Ingress_onset = tmp.Ingress_onset(randperm(640,640));
tmp.Tremble_Mag = tmp.Tremble_Mag(randperm(640,640));
tmp.Tremble_Dur = tmp.Tremble_Dur(randperm(640,640));
mdl = fitglm( tmp, modelspec, 'Distribution', 'normal' )

w = [0.261,-1.261,-0.1,0.005584]
predicted_ingress_mag = [ones(640,1),tmp.Tremble_Mag,tmp.Tremble_Dur,tmp.Ingress_onset]*w';

figure; scatter( tmp.Ingress_Mag, predicted_ingress_mag, 100, tmp.TrialNumber, 'filled' );
xlabel('Actual'); ylabel('Predicted')
xlim([0,1]);ylim([0,1]);line([0,1],[0,1])