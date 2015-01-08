%%
clear all
clf
colormap lines(11);
map = colormap;
leg = [];
ind = 1:2;
for i = ind
   name = ['Result_scripts\Commutativity\horse1-correspondancesAtError' num2str(i) '.mat']
   nInfo = ['Result_scripts\Commutativity\horse1-info' num2str(i) '.mat']
   a = load(name);
   a = a.correspondancesAtError;
   b = load(nInfo);
   b = b.infos;
   l{i-ind(1)+1} = b;
   hold on
   plot(linspace(0,40/383,4001),a/19248,'Color',map(i,:),'LineWidth',3);
end


legend(l)