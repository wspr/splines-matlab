function hobbysplines(points,varargin)
%HOBBYSPLINES Draws open or closed smooth curves from keypoints
%
% hobbysplines({[x1 y1], [x2 y2],... ,[xN,yN]},[opts])
%
% Draws a closed (cyclic) curve through keypoints 1 to N.
% Keypoints may be specified with optional slopes and tension parameters.
% The syntax to do this (replacing `[x1 y1]` above, say) is
%
%    { [x1 y1] s1 t1 }
%    { [x1 y1] s1 t1_in t1_out }
%
% where `s1` is the slope in degrees of the curve through `[x1 y1]`,
% and `t1_out` & `t1_in` are the "exit" and "entry" tensions of the curve
% approaching that point. If `t1` only is specified, this is both the
% entry and the exit tension.
%
% Use '' to indicate default values here; for example, for a default slope
% but specified "exit" tension of 2.0, use
%
%    { [x1 y1] '' 2.0 ''}
%
% and note that trailing optional arguments can also be omitted.
%
% Optional arguments given by [opts] can be any combination of the
% following:
%
%   OPTION       DEFAULT            DESCRIPTION
%   ------       -------            -----------
%   'tension'    [1]                default tension between points
%   'offset'     [0 0]              offset to add to each control point
%   'cycle'      [true]             draw an open curve instead
%   'debug'      [false]            draw and label keypoints on the curve
%   'linestyle'  '-'                line style option(s) (use a cell for multiple)
%   'color'      [0 0 0]            colour of the curve
%
% Note that at tension of 1 creates roughly circular plots

%% Parse inputs

p = inputParser;
p.addRequired('points',@iscell)
p.addOptional('tension',1);
p.addOptional('offset',[0 0]);
p.addOptional('cycle',true);
p.addOptional('color',[0 0 0]);
p.addOptional('debug',false);
p.addOptional('linestyle',{'color','black','linewidth',1});

p.parse(points,varargin{:});

cycle = p.Results.cycle;
offset = p.Results.offset;
debug = p.Results.debug;
points = p.Results.points;

color = p.Results.color;
linestyle = p.Results.linestyle;
if ~iscell(linestyle)
  linestyle = {linestyle};
end
  
if cycle
  points{end+1} = points{1};
end

Npoints = numel(points);

z = cell(Npoints,1);
w = cell(Npoints,1);
t = cell(Npoints,1);

for n = 1:Npoints
  
  pp = points{n};
  
  w{n} = [NaN NaN];
  tout{n} = p.Results.tension;
  tin{n} = p.Results.tension;
  
  if iscell(pp)
    
    veclen = numel(pp);
    z{n} = offset+pp{1};
    
    if veclen >= 2 && isnumeric(pp{2}) % lazy evaluation is my friend
      w{n} = [cosd(pp{2}) sind(pp{2})];
    end
    
    if veclen >= 3 && isnumeric(pp{3})
      tin{n} = pp{3};
    end
    
    if veclen == 4 && isnumeric(pp{4})
      tout{n} = pp{4};
    end
    
  else
    
    z{n} = offset+pp;
    
  end
  
end


%% fixup vectors iff necessary

if all( isnan(w{1}) )
  if cycle
    w{1} = z{2}-z{end-1};
  else
    w{1} = z{2}-z{1};
  end
  w{1} = w{1}/norm(w{1});
end
if all( isnan(w{end}) )
  if cycle
    w{end} = z{2}-z{end-1};
  else
    w{end} = z{end}-z{end-1};
  end
  w{end} = w{end}/norm(w{end});
end
for ii = 2:Npoints-1
  if all( isnan(w{ii}) )
    w{ii} = -z{ii-1} + z{ii+1};
  end
  w{ii} = w{ii}/norm(w{ii});
end

%% Calculate control points and plot bezier curve segments

hold on

for ii = 1:Npoints-1
  
  theta = arg(w{ii})-arg(z{ii+1}-z{ii});
  phi   = arg(z{ii+1}-z{ii})-arg(w{ii+1});
  
  [rho,sigma] = velocity_parameters(theta,phi);
  
  plot_bezier(...
    z{ii},...
    z{ii}+rho/(3*tout{ii})*norm(z{ii+1}-z{ii})*w{ii},...
    z{ii+1}-sigma/(3*tin{ii+1})*norm(z{ii+1}-z{ii})*w{ii+1},...
    z{ii+1},...
    {linestyle{:},'color',color})

end

if debug
  if cycle
    Mpoints = Npoints-1;
  else
    Mpoints = Npoints;
  end
  parfor ii = 1:Mpoints
    plot(z{ii}(1),z{ii}(2),'o','color',color)
    text(z{ii}(1),z{ii}(2),['   ',num2str(ii)])
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

function plot_bezier(P1,P2,P3,P4,linestyle)

N = 50;
t = linspace(0,1,N)';

c1 = 3*(P2 - P1);
c2 = 3*(P1 - 2*P2 + P3);
c3 = -P1 + 3*(P2-P3) + P4;

Q = t.^3*c3 + t.^2*c2 + t*c1 + repmat(P1,[N 1]);

plot(Q(:,1),Q(:,2),linestyle{:});

end

