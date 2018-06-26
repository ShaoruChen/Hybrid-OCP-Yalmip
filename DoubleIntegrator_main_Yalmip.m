clear;
T = 15;          % time horizon
d = 4;          % degree of relaxation
nmodes = 2;     % number of modes
r2 = 0.3;       % r^2, where r is the radius of the domain of mode 1

% Define variables
% t = msspoly( 't', 1 );
t = sdpvar(1);
x = cell( nmodes, 1 );
u = cell( nmodes, 1 );
f = cell( nmodes, 1 );
g = cell( nmodes, 1 );
x0 = cell( nmodes, 1 );
%%%%%%%%%%%%%%%%%%%%%%%%
% Can hX be expressed as a set of inequalities, rather than only one?
hX = cell( nmodes, 1 );
hU = cell( nmodes, 1 );
hXT = cell( nmodes, 1 );
sX = cell( nmodes, nmodes );
R = cell( nmodes, nmodes );
h = cell( nmodes, 1 );
H = cell( nmodes, 1 );

x0{2} = [ 1; 1 ];

% Dynamics
% sdpvar x1 x2;
% x{1} = [x1;x2];
x{1} = sdpvar(2,1,'full');
u{1} = sdpvar(1);
f{1} = T * [ x{1}(2); 0 ];
g{1} = T * [ 0; 1 ];

x{2} = x{1};
u{2} = u{1};
f{2} = f{1};
g{2} = g{1};

% Domains
% Mode 1 
y = x{1};
hX{1} = r2 - y'*y;
hU{1} = 1 - u{1}^2;
hXT{1} = hX{1};
h{1} = 1*x{1}'* x{1} + 20 * u{1}^2;
H{1} = 0;

% Mode 2
y = x{2};
hX{2} = y'*y - r2;
hU{2} = 1 - u{2}^2;
hXT{2} = hX{2};
sX{2,1} = [ r2 - y' * y;
            y' * y - r2 ];
R{2,1} = x{1};
h{2} = 1*x{2}' * x{2} + 20 * u{2}^2;
H{2} = 0;

% Options
options.MinimumTime = 0;
options.withInputs = 1;
options.svd_eps = 1e4;

[out] = HybridOCPDualSolver_Yalmip(t,x,u,f,g,hX,hU,sX,R,x0,hXT,h,H,d,options)
