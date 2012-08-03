%%

hfig = figure(1);
set(hfig,'color',[1 1 1],'name','Spline example');

clf;

points = {...
 [0 0],[1 1],[2,0],...
};
hobbysplines(points,'debug',true,'cycle',true);

points = {[1 0],[1 1],[0 1],[0 0]};
hobbysplines(points,'debug',true,'cycle',true);

axis equal
axis off

%%

hfig = figure(2); clf;
set(hfig,'color',[1 1 1],'name','Spline example');

points = {{[0 0] '' 1 1},[0.7 0.8],[0.8 2],[1 4],[0 5],[-1 3.5],[-0.8 2],[-0.8 1]};
NN = numel(points);
hobbysplines(points,'debug',true,'defaultTension',2,'cycle',true);

axis equal
axis off

