idx = 2;
Y = mymice.mouseObjs.PO340.displacement_csplus(idx,:);
ing = mymice.mouseObjs.PO340.ingress_index_csplus(idx);
Nframes = numel(Y);

figure; 
axes_ = arrayfun( @(x) subplot(3,1,x,'NextPlot','add'), [1:3] );
plot( axes_(1), Y );
d = 3000;
f0 = smooth(Y(ing:ing+d),5); plot( axes_(2), f0 );
f1 = -smooth(fliplr(Y(ing-d:ing)),5); plot( axes_(3), f1 );
[pks_max,locs_max] = findpeaks( f0, 'npeaks',1);
[pks_min,locs_min] = findpeaks( f1,'npeaks',1);
pks_min = -pks_min;

myline = linspace( pks_min, pks_max, numel( [ing - locs_min: ing + locs_max ] ) );

plot( axes_(1), ing-locs_min, pks_min, 'ro' );
plot( axes_(1), ing+locs_max, pks_max, 'ro' );

plot( axes_(2), locs_max, pks_max, 'ro' ); camroll(180)
plot( axes_(3), locs_min, -pks_min, 'ro' )
%plot( [ing - locs_min: ing + locs_max ], myline, 'r' )

% Before ingress fit
pc(1,:) = polyfit( [1:ing-locs_min] , Y([1:ing-locs_min]), 1);
% Ingress fit
pc(2,:) = [ (pks_max-pks_min)/(locs_max+locs_min) pks_min ];
% Post-ingress fit
pc(3,:) = polyfit( [ing+locs_max: Nframes] , Y([ing+locs_max:end]), 1);

plot( axes_(1), [1:ing-locs_min], polyval( pc(1,:), [1:ing-locs_min] ), 'k--' )
plot( axes_(1), [ing-locs_min:ing+locs_max], polyval( pc(2,:), [0:locs_min+locs_max] ), 'k--' )
plot( axes_(1), [ing+locs_max:Nframes], polyval( pc(3,:), [0:(Nframes-ing-locs_max)] ), 'k--' )

Fit_Table = table();
Fit_Table.X_pre = {1:ing-locs_min};
Fit_Table.Y_pre = {polyval( pc(1,:), [1:ing-locs_min] )};
Fit_Table.X_ing = {[ing-locs_min:ing+locs_max-1]};
Fit_Table.Y_ing = {polyval( pc(2,:), [0:locs_min+locs_max-1] )};
Fit_Table.X_post = {[ing+locs_max:Nframes-1]};
Fit_Table.Y_post = {polyval( pc(3,:), [0:(Nframes-ing-locs_max-1)] )};

%%

%%
%figure; cwt(Y(1:14000), 1000)
%

[WT,F] = cwt(Y,1000);
figure; imagesc(abs(WT))
%%
xrec = icwt( WT, F, [1 4] );
figure; plot( xrec )

%%

