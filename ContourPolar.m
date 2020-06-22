function hpol = ContourPolar(theta,rho,z,n,maxrho, SZA)
%ContourPolar Polar coordinate plot.
%   modified from POLAR in 2001/01/01
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.17 $  $Date: 1998/09/30 15:25:05 $
%theta 
%z 反射率值
%n 等值线的条数
%maho 极径的最大值


if nargin < 4
   n = 8;
end
if nargin < 5
   maxrho = max(rho);
end
if isstr(theta) | isstr(rho)
    error('Input arguments must be numeric.');
end
if ~isequal(size(theta),size(rho))
    error('THETA and RHO must be the same size.');
end

% get hold state
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

% transform data to Cartesian coordinates.
delta = maxrho/50; 
xx = rho.*cos((90-theta)*(pi/180));
yy = rho.*sin((90-theta)*(pi/180));
xi = -maxrho:delta:maxrho;
yi = xi;
[xi,yi] = meshgrid(xi,yi);
zi = griddata(xx,yy,z,xi,yi,'cubic');
contourf(xi,yi,zi,n);
%surf(xi,yi,zi);
%axis([-60,60,-60,60,25,40]);
%view(0,90);
%shading flat;


% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');
ls = get(cax,'gridlinestyle');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
% fSize = 18;
% fName = 'Time New Roman';
fWeight = get(cax, 'DefaultTextFontWeight');
fUnits  = get(cax, 'DefaultTextUnits');
set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
    'DefaultTextFontName',   get(cax, 'FontName'), ...
    'DefaultTextFontSize',   get(cax, 'FontSize'), ...
    'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
    'DefaultTextUnits','data')

% only do grids if hold is off
if ~hold_state
% make a radial grid
	 hold on;
%    maxrho = max(abs(rho(:)));
    hhh=plot([-maxrho -maxrho maxrho maxrho],[-maxrho maxrho maxrho -maxrho]);
    set(gca,'dataaspectratio',[1 1 1],'plotboxaspectratiomode','auto')
    v = [get(cax,'xlim') get(cax,'ylim')];
    ticks = sum(get(cax,'ytick')>=0);
    delete(hhh);
% check radial limits and ticks
    rmin = 0; rmax = v(4); rticks = max(ticks-1,2);
    if rticks > 5   % see if we can reduce the number
        if rem(rticks,2) == 0
            rticks = rticks/2;
        elseif rem(rticks,3) == 0
            rticks = rticks/3;
        end
    end

% define a circle
    th = 0:pi/50:2*pi;
    xunit = cos(th);
    yunit = sin(th);
% now really force points on x/y axes to lie on them exactly
    inds = 1:(length(th)-1)/4:length(th);
    xunit(inds(2:2:4)) = zeros(2,1);
    yunit(inds(1:2:5)) = zeros(3,1);
% plot background if necessary
%    if ~isstr(get(cax,'color')),
%       patch('xdata',xunit*rmax,'ydata',yunit*rmax, ...
%             'edgecolor',tc,'facecolor',get(gca,'color'),...
%             'handlevisibility','off');
%    end

% draw radial circles
    c82 = cos(82*pi/180);
    s82 = sin(82*pi/180);
    rinc = (rmax-rmin)/rticks;
    for i=(rmin+rinc):rinc:rmax
        hhh = plot(xunit*i,yunit*i,ls,'color',tc,'linewidth',1,...
                   'handlevisibility','off');
        text((i+rinc/20)*c82,(i+rinc/20)*s82, ...
            ['  ' num2str(i)],'verticalalignment','bottom',...
            'handlevisibility','off','fontname','Time New Roman','fontsize',6)
    end
    set(hhh,'linestyle','-') % Make outer circle solid

% plot spokes
    th = (1:2)*2*pi/4;
    cst = cos(pi/2-th); snt = sin(pi/2-th);
    cs = [-cst; cst];
    sn = [-snt; snt];
    plot(rmax*cs,rmax*sn,ls,'color',tc,'linewidth',0.8,...
         'handlevisibility','off')

% annotate spokes in degrees
    rt = 1.15*rmax;
    for i = 1:length(th)
        text(rt*cst(i),rt*snt(i),int2str(i*90),...
             'horizontalalignment','center',...
             'handlevisibility','off','fontname','Time New Roman','fontsize',6);
        if i == length(th)
            loc = int2str(0);
        else
            loc = int2str(180+i*90);
        end
        text(-rt*cst(i),-rt*snt(i),loc,'horizontalalignment','center',...
             'handlevisibility','off','fontname','Time New Roman','fontsize',6)
    end

% set view to 2-D
    view(2);
% set axis limits
    axis(rmax*[-1 1 -1.15 1.15]);
end

% Reset defaults. fSize
set(cax, 'DefaultTextFontAngle', fAngle , ...
    'DefaultTextFontName',   fName , ...
    'DefaultTextFontSize',   fSize , ...
    'DefaultTextFontWeight', fWeight, ...
    'DefaultTextUnits',fUnits );
% plot data on top of grid
if nargout > 0
    hpol = q;
end
if ~hold_state
    set(gca,'dataaspectratio',[1 1 1]), axis off; set(cax,'NextPlot',next);
end

hold on
plot(xx,yy,'.k', 'MarkerSize', 10);
plot(0, SZA, 'pr', 'MarkerSize', 6, 'MarkerFaceColor', 'r');

set(get(gca,'xlabel'),'visible','on');
set(get(gca,'ylabel'),'visible','on');
colormap gray;
%pcolor(255,255,255);


