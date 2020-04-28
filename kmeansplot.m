ingress_index_csplus
ingress_index_csminus

% Need max of displacement_csplus and displacement_csminus

trembleScore_csplus
trembleScore_csminus

ingressScore_csplus
ingressScore_csminus

ingress_stats_CSplus
ingress_stats_CSminus [2:3]


%%

% CS+ Lat, CS+ Score, Ingress Fit coeff1, Ingress Fit coeff2, Tremble Score
forkmeans_csplus = [];
mouse = fields(mymice.mouseObjs)';
for i = 1:length(mouse)
   i0 = mymice.mouseObjs.(mouse{i}).ingress_index_csplus';
   i1 = mymice.mouseObjs.(mouse{i}).ingressScore_csplus';
   i2 = cellfun(@(x) x(2), mymice.mouseObjs.(mouse{i}).ingress_stats_CSplus, 'ErrorHandler', @(a,b) NaN);
   i3 = cellfun(@(x) x(3), mymice.mouseObjs.(mouse{i}).ingress_stats_CSplus, 'ErrorHandler', @(a,b) NaN);
   i4 = mymice.mouseObjs.(mouse{i}).trembleScore_csplus';

   forkmeans_csplus = [forkmeans_csplus;[i0 i1 i2 i3 log(i4)]];

end
    forkmeans_csplus(:,6)=reshape(repmat([1:19],10,1),190,1);
    forkmeans_csplus(:,7)=repmat([1:10],1,19)';
%%
forkmeans_csminus = [];
for i = 1:length(mouse)
   i0 = mymice.mouseObjs.(mouse{i}).ingress_index_csminus';
   i1 = mymice.mouseObjs.(mouse{i}).ingressScore_csminus';
   i2 = cellfun(@(x) x(2), mymice.mouseObjs.(mouse{i}).ingress_stats_CSminus, 'ErrorHandler', @(a,b) NaN);
   i3 = cellfun(@(x) x(3), mymice.mouseObjs.(mouse{i}).ingress_stats_CSminus, 'ErrorHandler', @(a,b) NaN);
   i4 = mymice.mouseObjs.(mouse{i}).trembleScore_csminus';
   
   forkmeans_csminus = [forkmeans_csminus;[i0 i1 i2 i3 i4]];
   
end


%%
nanz = @(x) (x-nanmean(x))/nanstd(x);
for i=1:5; 
   
    z_forkmeans_csplus(:,i)=nanz(forkmeans_csplus(:,i));
    
end

%%
[idx,ctrs]=kmeans(z_forkmeans_csplus,2);

%%
figure;
dm=[2,5]; %ingress score and tremble score
plot(z_forkmeans_csplus(idx==1,dm(1)),z_forkmeans_csplus(idx==1,dm(2)),'r.','MarkerSize',12);
hold on;
plot(z_forkmeans_csplus(idx==2,dm(1)),z_forkmeans_csplus(idx==2,dm(2)),'b.','MarkerSize',12);
plot(ctrs(:,2),ctrs(:,5),'ko','MarkerSize',12,'LineWidth',2);
legend('Cluster 1','Cluster 2','Centroids');
hold off;
%%
img = mymice.mouseObjs.PO341.cwt_csplus{5};
img = real(img([200:350],2000:6000));
[Gx,Gy]= imgradientxy( img );
min_Gy = graythresh(Gy(:));
x = arrayfun( @(x) (x<min_Gy), Gy );
max(max( x.*img))


%%

figure; subplot(4,1,1); imagesc(img);
subplot(4,1,2); imagesc(grd);
subplot(4,1,3); imagesc(x);
subplot(4,1,4); imagesc( x.*img );

%%


figure;

colors = lines(19);

ax = arrayfun(@(x) subplot(5,4,x, 'xlim', [0,1], 'ylim', [-7,0], 'NextPlot','add'), [1:20]); 
shape = 'square'
arrayfun( @(animal) plot( forkmeans_csplus( find(forkmeans_csplus(:,6)==animal), 2),...
    forkmeans_csplus(find(forkmeans_csplus(:,6)==animal), 5 ) , shapes{animal},...
    'MarkerFaceColor', colors(animal,:), 'MarkerSize', 12, 'Parent', ax(animal) ), [1:19] )

mousenames = fields(mymice.mouseObjs);
arrayfun( @(animal) text(0,0,sprintf('%s',mousenames{animal}),'parent',ax(animal)), [1:19])

%%

trembleduration = structfun( @(x) arrayfun( @(y) sum( max( x.cwt_csplus{y}([250:350],[2000: x.ingress_index_csplus(y) ]),[],1 ) >(5*10^-4) ), [1:10], 'uniformoutput',false ), mymice.mouseObjs, 'uniformoutput',false );

%%

my_duration_fields = fields(trembleduration);

tremble_duration_column = [];
for i = 1:numel( my_duration_fields )
   
    tremble_duration_column = [tremble_duration_column; ...
        (arrayfun( @(x) trembleduration.(my_duration_fields{i}){x}, [1:10] ))' ];
    
end

forkmeans_csplus(:,8) = tremble_duration_column;

%% 


figure;

colors = lines(19);

ax = arrayfun(@(x) subplot(5,4,x, 'xlim', [-7,0], 'ylim', [0,12000], 'NextPlot','add'), [1:20]); 
shape = 'square'
arrayfun( @(animal) plot( forkmeans_csplus( find(forkmeans_csplus(:,6)==animal), 5),...
    forkmeans_csplus(find(forkmeans_csplus(:,6)==animal), 8 ) , shapes{animal},...
    'MarkerFaceColor', colors(animal,:), 'MarkerSize', 12, 'Parent', ax(animal) ), [1:19] )

mousenames = fields(mymice.mouseObjs);
arrayfun( @(animal) text(0,0,sprintf('%s',mousenames{animal}),'parent',ax(animal)), [1:19])
