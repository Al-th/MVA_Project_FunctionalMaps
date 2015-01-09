%%
clear all
clf
colormap lines(11);
map = colormap;
leg = [];
ind = 1:10;
figure(1);
hold on
for i = 1:2
   name = ['Result_scripts\2_Refinement\geo_error.0003.isometry.5_' num2str(i-1) '.mat']
   a = load(name);
   a = a.pts;
   l{i} = ['Refinement ' num2str(i-1)];
   plot(linspace(0,40/383,4001),a/19248,'Color',map(i,:),'LineWidth',3);
end
xlabel('Normalized geodesic error') % x-axis label
ylabel('Percentage of true correspondences') % y-axis label
legend(l)