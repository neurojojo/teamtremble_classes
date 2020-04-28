%Download https://drive.google.com/open?id=1rB26y0EkQ7c-LFinIVVAVvj_Cgd06xcf
%Place ActivePassiveAnalysis folder in 'C:\ActivePassiveAnalysis\'
clear all

addpath( genpath('C:\teamtremble_classes') )

mymice = findDisplacedMice( 'C:\ActivePassiveAnalysis' );

%
mymice.associateKey('C:\\teamtremble_classes\\cs_all_indices_key.csv');

mymice.loadDisplacements();

mymice.makeList()
% Manually curated file with three columns (Animal ID, third index of
% displacement for CS+, third index of displacement for CS-)


mymice.associateIngressTimes('C:\\teamtremble_classes\\ingress_times_noingress_set_to_10000.csv');

mymice.getIngressScore_allMice();
mymice.fitLines_allMice('ingress');
mymice.fitLines_allMice('trial');
mymice.getTrembleScore_allMice();


