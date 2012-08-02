function hobbysplines(points,varargin)
%HOBBYSPLINES Draws open or closed smooth curves from keypoints
%
% hobbysplines({[x1 y1], [x2 y2],... ,[xN,yN]},[opts])
%
% Draws a closed (cyclic) curve through keypoints 1 to N.
% Use the ['cycle',false] option to draw an open curve instead.

%% Parse inputs

p = inputParser;
p.addRequired('points',@iscell)
p.addOptional('slopes',false);
p.addOptional('tensions',false);
p.addOptional('cycle',true);
p.addOptional('debug',false);
p.addOptional('linestyle',{'color','black','linewidth',1});

p.parse(points,varargin{:});

z = p.Results.points;
cycle = p.Results.cycle;
debug = p.Results.debug;

if numel(p.Results.slopes) > 1
  w = arrayfun(@(x) [cosd(x) sind(x)],p.Results.slopes,'UniformOutput',false);
else
  w = repmat({[0 0]},size(z));
end

if iscell(p.Results.tensions)
  t = p.Results.tensions;
else
  t = repmat({[0 0]},size(z));
end


if cycle
  z{end+1} = z{1};
  if numel(w) ~= numel(z)
    w{end+1} = w{1};
  end
  if numel(t) ~= numel(z)
    t{end+1} = [1 1];
  end
end

Npoints = numel(z);

%% fixup vectors iff necessary

t(cellfun(@(x) all(x==[0 0]),t)) = {[1 1]};

if all( w{1}==[0 0] )
  if cycle
    w{1} = z{2}-z{end};
  else
    w{1} = z{2}-z{1};
  end
end
if all( w{end}==[0 0] )
  if cycle
    w{end} = z{2}-z{end};
  else
    w{end} = z{end}-z{end-1};
  end
end
for ii = 2:Npoints-1
  if all( w{ii}==[0 0] )
    w{ii} = -z{ii-1} + z{ii+1};
  end
end

%% Calculate control points and plot bezier curve segments
%
% My bezier code (stolen, below) can't accept more than four control points
% and still pass through the desired points of each. Not sure whether this
% is intentional.

hold on

for ii = 1:Npoints-1
  
  theta = arg(w{ii})-arg(z{ii+1}-z{ii});
  phi   = arg(z{ii+1}-z{ii})-arg(w{ii+1});
  
  [rho,sigma] = velocity_parameters(theta,phi);
  
  control_points = [...
    z{ii};
    z{ii}+rho/(3*t{ii}(2))*norm(z{ii+1}-z{ii})*w{ii};
    z{ii+1}-sigma/(3*t{ii+1}(1))*norm(z{ii+1}-z{ii})*w{ii+1};
    z{ii+1};
    ];
  
  mybez(control_points,p.Results.linestyle)

end

if debug
  parfor ii = 1:Npoints
    plot(z{ii}(1),z{ii}(2),'o','color',[1 0 0])
  end
end

end

%% Sub-functions

function o = arg(w)
  o = atan2(w(2),w(1));
end

function [rho,sigma] = velocity_parameters(theta,phi)
% From "Smooth, easy to compute interpolating splines" by John D. Hobby
% <http://www.springerlink.com/content/p4w1k8w738035w80/>

a = 1.597;
b = 0.07;
c = 0.37;

st = sin(theta);
ct = cos(theta);
sp = sin(phi);
cp = cos(phi);

alpha = a*(st-b*sp)*(sp-b*st)*(ct-cp);
rho   = (2+alpha)/(1+(1-c)*ct+c*cp);
sigma = (2-alpha)/(1+(1-c)*cp+c*ct);

end

function mybez(P,linestyle)
% This function is an edited version of mybez.m by Sagar Aiya.
% Original source, licensed under BSD:
% <http://www.mathworks.com/matlabcentral/fileexchange/30759-bezier-curve-plotter/>
%
% This is probably overkill for what I need but it was the first thing I
% stumbled up :) I'd prefer to reimplement something myself if time
% permits.

P = transpose(P);
n = length(P);
count = 1;

div = 50; %number of segments of the curve (Increase this value to obtain a
          %smoother curve

for u = 0:(1/div):1
    sum = [0 0]';
    for i = 1:n
        B = nchoosek(n,i-1)*(u^(i-1))*((1-u)^(n-i+1)); %B is the Bernstein polynomial value
        sum = sum + B*P(:,i);
    end
    B = nchoosek(n,n)*(u^(n));
    sum = sum + B*P(:,n);
    A(:,count) = sum; %the matrix containing the points of curve as column vectors. 
    count = count+1;  % count is the index of the points on the curve.
end

x = A(1,:);
y = A(2,:);
plot(x,y,linestyle{:});

end

