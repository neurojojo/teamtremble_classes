%Download https://drive.google.com/open?id=1rB26y0EkQ7c-LFinIVVAVvj_Cgd06xcf
%Place ActivePassiveAnalysis folder in 'C:\teamtremble\data\'

mymice = findDisplacedMice( 'C:\teamtremble\data\ActivePassiveAnalysis' );

mymice.associateKey('C:\\teamtremble\\cs_third_index_key.csv');
mymice.associateIngressTimes('C:\\teamtremble\\ingress_times.csv');

mymice.loadDisplacements();

mymice.computeFourier();