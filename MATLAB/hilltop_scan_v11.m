
%% hilltop_scan_v11
% 13 Jan 2025
% Pete and Claire discussions
%
clear all;
%close all;

%% New params post SIADS discussion 15 May 2025
% load parameter array
hilltop_scan_params6

figure;
clf;
f=gcf();
f.Position(3:4)=[950 450];

cmap=colormap('parula');
cmap2=flipud(cmap);

tic;
% edit this if you want to select which are generated
%for ii=7:nps
for ii=1:nps
%for ii=4
    clf;
    p=ps(ii);
    
    npts=150;
    [X,Y]=meshgrid(p.smin:(p.smax-p.smin)/npts:p.smax,p.smin:(p.smax-p.smin)/npts:p.smax);

    % x nullcline
    xa=X(1,:);
    ya=nullfory(xa,p);

    % y nullcline
    yb=Y(:,1)';
    xb=nullforx(yb,p);


    %% compute escape time and angle
    p.tmax=40;
    TT=[0 p.tmax];
    ZEND=zeros(size(X));
    TEND=zeros(size(X));
    odeopt = odeset('Events',@(t,y) leavebox(t,y,p));
    for i=1:npts+1
        for j=1:npts+1
            xx=[X(i,j);Y(i,j)];
            [tt,yy]=ode45(@(t,x) hilltop_fn(x,p),TT,xx,odeopt);
            TEND(i,j)=tt(end);
            ZEND(i,j)=atan2(yy(end,2)-p.smin,yy(end,1)-p.smin)*4/pi-1;
            % test for no escape
            if yy(end,1)<p.smax && yy(end,2)<p.smax
                TEND(i,j)=NaN;
                ZEND(i,j)=NaN;
            end
        end
    end

    % escape angle
    subplot(1,2,1)
    s1=pcolor(X,Y,ZEND);
    axis square
    hold on
    s1.EdgeColor = 'none';
    colourmap=turbo(20);
    colormap(gca,colourmap)
    %lcm=lines(200);
    %colormap(gca,lcm);
    clim([-1 1]);
    plot(xa,ya,'k');
    plot(xb,yb,'k');
    box on
    set(gca, 'Layer', 'Top');

    % escape time
    subplot(1,2,2)
    s2=pcolor(X,Y,real(TEND));
    s2.EdgeColor = 'none';
    axis square
    colormap(gca,cmap2);
    %colorbar
    clim([0 5]);
    hold on
    % nullclines
    plot(xa,ya,'k');
    plot(xb,yb,'k');
    box on
    set(gca, 'Layer', 'Top');
    drawnow();
    
  %     save for parameter set ii
    fname=sprintf('fig_abp-bottle-mu-%f-gamma-%f-%s.pdf',p.mu,p.gamma,p.name);
  %  toc;
    exportgraphics(gcf, fname, 'ContentType', 'vector');
end

%% plot only colorbars
%
% 
% figure(2);
% clf;
% f=gcf();
% f.Position(3:4)=[400 150];
% clim([0 100])
% cmap=colormap('parula');
% c = colorbar("north");
% set(gca,'Visible',false)
% c.Position = [0.1 0.3 0.74 0.5];
% 
% exportgraphics(f, "fig-colorbar1.pdf", 'ContentType', 'vector');
% 
% %%
% 
% figure(3);
% clf;
% f=gcf();
% f.Position(3:4)=[400 150];
% clim([-1 1])
% colormap(centered('RdYlGn'));
% c = colorbar("north");
% set(gca,'Visible',false)
% c.Position = [0.1 0.3 0.74 0.5];
% 
% exportgraphics(f, "fig-colorbar2.pdf", 'ContentType', 'vector');
% 
% %%
% 
% figure(4);
% clf;
% f=gcf();
% f.Position(3:4)=[400 150];
% clim([0 5])
% %colormap(flipud('parula'));
% cmap=colormap('parula');
% colormap(flipud(cmap));
% c = colorbar("north");
% set(gca,'Visible',false)
% c.Position = [0.1 0.3 0.74 0.5];
% 
% exportgraphics(f, "fig-colorbar3.pdf", 'ContentType', 'vector');
% 
% 
% %%


keyboard
%% function definitions

function [value, isterminal, direction]=leavebox(t,y,p)
value=max([y(1)-p.smax,y(2)-p.smax]);
isterminal=1;
direction=0;
end


function zz=hilltop_fn(xx,p)
x=xx(1,:);
y=xx(2,:);
alpha=p.alpha;
beta=p.beta;
gamma=p.gamma;
mu=p.mu;
zz=zeros(size(xx));
zz(1,:)=x.^2-(mu+gamma*(gamma+alpha^2-beta^2)/((alpha+beta)^2))+2*alpha*y-gamma;
zz(2,:)=y.^2-(mu+gamma*(gamma+alpha^2-beta^2)/((alpha+beta)^2))+2*beta*x+gamma;
end

function zz=hilltop_jac(xx,p)
x=xx(1);
y=xx(2);
alpha=p.alpha;
beta=p.beta;
zz=zeros(2,2);
zz(1,1)=2*x;
zz(1,2)=2*alpha;
zz(2,1)=2*beta;
zz(2,2)=2*y;
end

function y=nullfory(x,p)
alpha=p.alpha;
beta=p.beta;
gamma=p.gamma;
mu=p.mu;
y=-(x.^2-(mu+gamma*(gamma+alpha^2-beta^2)/((alpha+beta)^2))-gamma)/(2*alpha);
end

function x=nullforx(y,p)
alpha=p.alpha;
beta=p.beta;
gamma=p.gamma;
mu=p.mu;
x=-(y.^2-(mu+gamma*(gamma+alpha^2-beta^2)/((alpha+beta)^2))+gamma)/(2*beta);
end


