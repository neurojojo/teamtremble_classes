mouse = 'PO374';
trial = 5;

figure('color','w'); imagesc( log(mymice.mouseObjs.(mouse).cwt_csplus{trial}(:,:)) ); colormap( hot(300) )
set(gca,'Clim',[-11,-2],'TickDir','out','XTick',[],'YTick',[])
hold on; plot( 100-50*smooth(mymice.mouseObjs.(mouse).displacement_csplus(trial,:),5), 'c', 'LineWidth', 2 ); 

ingressBuffer=1000

line([2000,mymice.mouseObjs.(mouse).ingress_index_csplus(trial)-ingressBuffer],[5,5],'color','w','linewidth',4)
line([mymice.mouseObjs.(mouse).ingress_index_csplus(trial)-ingressBuffer,mymice.mouseObjs.(mouse).ingress_index_csplus(trial)], [180,180], 'color', 'w','linewidth',4)

saveas(gcf, sprintf('representative_%s_%i.svg', mouse, trial ) )
%%

tmp = mymice.mouseObjs.PW28;
i = 10
ingressBuffer = 500

for i = 1:20
result(i) = tmp.max_of_grad_adj_img( tmp.cwt_csplus{i}(:,[2000: tmp.ingress_index_csplus(i)-ingressBuffer ]))
end