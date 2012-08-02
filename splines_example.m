
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
