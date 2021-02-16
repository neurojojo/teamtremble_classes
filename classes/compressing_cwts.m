% To add to the demo script %
% January 5th, 2021 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%%

% Here we decrease the size of the CWT by converting it to uint8
flatten = @(x) x(:);
rescale = @(x,min,max) (x-min(flatten(x(:))))./(max(flatten(x(:)))-min(flatten(x(:))));
%converted_ = arrayfun( @(x) convertToUint8( x{1} ), mymice.mouseObjs.PO330.cwt_csplus, 'UniformOutput', false );
max_ = [];
for field = fields(mymice.mouseObjs)'
    max_ = [ max_ ; arrayfun( @(x) max(flatten(x{1})), mymice.mouseObjs.(field{1}).cwt_csplus ) ];
end

%% Compressed structure
mymice_not_object_cwt = struct();
mymice_not_object_displacement = struct();

for field = fields(mymice.mouseObjs)'
    mymice_not_object_cwt.(field{1}).cwt_csplus = arrayfun( @(x) convertToUint16(x{1},0,100), mymice.mouseObjs.(field{1}).cwt_csplus, 'UniformOutput', false );
    mymice_not_object_cwt.(field{1}).cwt_csminus = arrayfun( @(x) convertToUint16(x{1},0,100), mymice.mouseObjs.(field{1}).cwt_csminus, 'UniformOutput', false );
    mymice_not_object_displacement.(field{1}).displacement_csplus = arrayfun( @(x) convertToUint16(x,0,20), mymice.mouseObjs.(field{1}).displacement_csplus, 'UniformOutput', true );
    mymice_not_object_displacement.(field{1}).displacement_csminus = arrayfun( @(x) convertToUint16(x,0,20), mymice.mouseObjs.(field{1}).displacement_csminus, 'UniformOutput', true );
end

%%
tmp = convertToUint16( mymice.mouseObjs.PO330.cwt_csplus{1},0,100 );
%% Arranged as 10 trials csplus 1, 10 trials csminus 1, 10 trials csplus 2, 10 trials csminus 2
ingressTable = table2array( mymice.ingressTable(:,[2:end]) );
%% Saving files

save('mymice_not_object_displacement.mat','mymice_not_object_displacement')
save('mymice_not_object_cwt.mat','mymice_not_object_cwt')

ingressTable_labels = table2array( ingressTable.Labels(:,1) ); ingressTable_values = table2array(mymice.ingressTable(:,[2:end]));
save('mymice_ingresstable_labels.mat','ingressTable_labels')
save('mymice_ingressTable_values.mat','ingressTable_values')

%%
function output = convertToUint16( variable, min_, max_ )

    flatten = @(x) x(:);
    rescale = @(x,min_,max_) (x-min_)./(max_-min_);
    output = uint16( rescale( variable, min_, max_ ) * 2^16 );

end
