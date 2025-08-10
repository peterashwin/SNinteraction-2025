%% hilltop_scan_params
%%%%%%%%% Old->New renumbering 15th May 2025
% a a
% b b
% c c
% c’  d
% d f
% e g
% f h
% g i
% h j
% i k
% j l
% k e
%%%%%%%%%%

% default plotting parameters
p.smin=-10;
p.smax=10;

% default system parameters
p.alpha=3.1;
p.beta=-1.3;
p.mu=0.0;
p.gamma=0.0;
p.name="default";

pno=0;

%
pno=pno+1;
ps(pno)=p;
ps(pno).mu=-0.7;
ps(pno).gamma=4.1;
ps(pno).name="2a";


%
pno=pno+1;
ps(pno)=p;
ps(pno).mu=-0.7;
ps(pno).gamma=3.6;
ps(pno).name="2b";


%
pno=pno+1;
ps(pno)=p;
ps(pno).mu=-0.7;
ps(pno).gamma=-2.5;
ps(pno).name="2c";


%
pno=pno+1;
ps(pno)=p;
ps(pno).mu=-0.7;
ps(pno).gamma=-5;
ps(pno).name="2d";

%
pno=pno+1;
ps(pno)=p;
ps(pno).mu=-0.4;
ps(pno).gamma=1.65;
ps(pno).name="2e";

%
pno=pno+1;
ps(pno)=p;
ps(pno).mu=-0.2;
ps(pno).gamma=3.0;
ps(pno).name="2f";


%
pno=pno+1;
ps(pno)=p;
ps(pno).mu=-0.02;
ps(pno).gamma=0.69;
ps(pno).name="2g";

%
pno=pno+1;
ps(pno)=p;
ps(pno).mu=0.2;
ps(pno).gamma=3.5;
ps(pno).name="2h";

%
pno=pno+1;
ps(pno)=p;
ps(pno).mu=0.2;
ps(pno).gamma=-2.0;
ps(pno).name="2i";


%
pno=pno+1;
ps(pno)=p;
ps(pno).mu=0.4;
ps(pno).gamma=-2.0;
ps(pno).name="2j";

%
pno=pno+1;
ps(pno)=p;
ps(pno).mu=0.005;
ps(pno).gamma=3.58;
ps(pno).name="2k";

%
pno=pno+1;
ps(pno)=p;
ps(pno).mu=-0.005;
ps(pno).gamma=3.58;
ps(pno).name="2l";




%
pno=pno+1;
ps(pno)=p;
ps(pno).alpha=3.1;
ps(pno).beta=1.3;
ps(pno).mu=-10.0;
%ps(pno).mu=-6.0;
ps(pno).gamma=10.0;
ps(pno).name="1a";

%
pno=pno+1;
ps(pno)=p;
ps(pno).alpha=3.1;
ps(pno).beta=1.3;
ps(pno).mu=5.0;
ps(pno).gamma=10.0;
ps(pno).name="1b";


%
pno=pno+1;
ps(pno)=p;
ps(pno).alpha=3.1;
ps(pno).beta=1.3;
ps(pno).mu=17.0;
ps(pno).gamma=10.0;
ps(pno).name="1c";

%
pno=pno+1;
ps(pno)=p;
ps(pno).alpha=3.1;
ps(pno).beta=1.3;
ps(pno).mu=-10.0;
ps(pno).gamma=-13.0;
ps(pno).name="1d";

nps=pno;