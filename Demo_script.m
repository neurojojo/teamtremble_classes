%Download https://drive.google.com/open?id=1rB26y0EkQ7c-LFinIVVAVvj_Cgd06xcf
%Place ActivePassiveAnalysis folder in 'C:\ActivePassiveAnalysis\'

mymice = findDisplacedMice( 'C:\ActivePassiveAnalysis' );

%%
addpath( genpath('C:\teamtremble_classes') )


% Manually curated file with three columns (Animal ID, third index of
% displacement for CS+, third index of displacement for CS-)
mymice.associateKey('C:\\teamtremble_classes\\cs_third_index_key.csv');

mymice.associateIngressTimes('C:\\teamtremble_classes\\ingress_times_noingress_set_to_10000.csv');

mymice.loadDisplacements();

mymice.getTrembleScore_allMice();

mymice.getIngressScore_allMice();

%%


