%clear all

addpath( genpath('C:\teamtremble_classes') )
mymice_dreadd = findDisplacedMice( 'C:\ActivePassiveAnalysis2\DREADD inactivation\' );

mymice_dreadd.associateKey('C:\\teamtremble_classes\\cs_all_indices_key_dreadd.csv');

mymice_dreadd.loadDisplacements();


mymice_dreadd.makeList()

%

mymice_dreadd.associateIngressTimes('C:\\teamtremble_classes\\ingress_times_noingress_set_to_10000_dreadd.csv');

%

mymice_dreadd.fitLines_allMice('coeffs');
mymice_dreadd.fitLines_allMice('fitlines');
mymice_dreadd.getIngressScore_allMice();
mymice_dreadd.getTrembleScore_allMice();

%% Only do this once  %%%%%%
% 
% to_table = reshape( cell2mat( struct2cell ( structfun( @(x) x.ingressTable.Changepoint, mymice.mouseObjs , 'UniformOutput', false) ) ), 40, 6 )';
% 
% csplus = to_table(:,[1:20]);
% csminus = to_table(:,[21:end]);
% 
% colnames = [ {'matched_names'}, arrayfun( @(x) sprintf('ingress_csplus_matched_to_name_%i',x) , [1:20] , 'UniformOutput', false),...
% arrayfun( @(x) sprintf('ingress_csminus_matched_to_name_%i',x) , [1:20] , 'UniformOutput', false)]
% 
% all_ingress_table = [cell2table( fields(mymice.mouseObjs) ), array2table(csplus), array2table(csminus)];
% 
% all_ingress_table.Properties.VariableNames = colnames
% writetable( all_ingress_table, 'ingress_times_noingress_set_to_10000_dreadd.csv' )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%


mymice_dreadd.associateIngressTimes('C:\\teamtremble_classes\\ingress_times_noingress_set_to_10000_dreadd.csv');


mymice_dreadd.fitLines_allMice('coeffs');
mymice_dreadd.fitLines_allMice('fitlines');
mymice_dreadd.getIngressScore_allMice();
mymice_dreadd.getTrembleScore_allMice();

%%

