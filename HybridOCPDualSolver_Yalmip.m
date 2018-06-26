function [out] = HybridOCPDualSolver_Yalmip(t,x,u,f,g,hX,hU,sX,R,x0,hXT,h,H,d,options)


% Sanity check
if mod(d,2) ~= 0
    warning('d is not even. Using d+1 instead.');
    d = d+1;
end
nmodes = length(x);

max_m = 0;
for i = 1 : nmodes
    m = length(u{i});
    if m > max_m
        max_m = m;
    end
    if (length(f{i}) ~= length(x{i})) || (size(g{i},1) ~= length(x{i}))
        error('Inconsistent matrix size.');
    end
    if size(g{i},2) ~= m
        error('Inconsistent matrix size.');
    end
end

svd_eps = 1e3;
if nargin < 15, options = struct(); end
if isfield(options, 'svd_eps'), svd_eps = options.svd_eps; end
if ~isfield(options, 'freeFinalTime'), options.freeFinalTime = 0; end
if ~isfield(options, 'withInputs'), options.withInputs = 0; end

T = 1;
hT = t*(T-t);

if isempty(R)
    R = cell(nmodes,nmodes);
    for i = 1 : nmodes
    for j = 1 : nmodes
        if ~isempty(sX{i,j}), R{i,j} = x{i}; end
    end
    end
end

%% Setup spotless program
% define the program

% Create program variables in each mode
for i = 1 : nmodes
    
%     prog = prog.withIndeterminate( x{i} );
%     prog = prog.withIndeterminate( u{i} );
    
    % create v(i)
%   vmonom{ i } = monomials( [ t; x{ i } ], 0:d );
%   [ prog, v{ i }, ~ ] = prog.newFreePoly( vmonom{ i } );
    [v{i}, coef{i}, vmono{i}]  = polynomial([t; x{i}], d);
    % create the variables that will be used later
    vT{ i } = replace( v{ i }, t, T );
    dvdt{ i } = jacobian( v{ i }, t );
    dvdx{ i } = jacobian( v{ i }, x{ i } );
    Lfv{ i } = dvdt{ i } + dvdx{ i } * f{ i };
    Lgv{ i } = dvdx{ i } * g{ i };
    Lv{ i } = Lfv{ i } + Lgv{ i } * u{ i };
end

% creating the constraints and cost function
obj = 0;
F = [];  c = [];
for i = 1:nmodes
    c = [c; coef{i}];
end

for i = 1 : nmodes
    % Lv_i + h_i >= 0                   Dual: mu

    [F,c] = sosOnK_Yalmip(F,Lv{ i } + h{ i }, [ t; x{ i }; u{ i } ], [ hT; hX{ i }; hU{ i } ], d, c);
   
%     mu_idx(i) = size( prog.sosExpr, 1 );
    
    % v(T,x) <= H_i(x)                  Dual: muT
    if ~isempty( hXT{ i } )
        if options.freeFinalTime
            [F,c] = sosOnK_Yalmip( F, H{ i } - v{ i }, [ t; x{ i } ], [ hT; hXT{ i } ], d , c);
        else
            [F,c] = sosOnK_Yalmip( F, H{ i } - vT{ i }, x{ i }, hXT{ i }, d , c);
        end
    end
    
    % v_i(t,x) <= v_j (t,R(x))          Dual: muS
    for j = 1 : nmodes
        if ( ~isempty( sX{ i, j } ) ) % if its empty there isn't a guard between these
            vj_helper = replace( v{ j }, x{ j }, R{ i, j } ); 
            [F,c] = sosOnK_Yalmip( F, vj_helper - v{ i }, [ t; x{ i } ] , [ hT; sX{ i, j } ], d , c);
        end
    end
    
    % Objective function
    if ~isempty( x0{ i } )
        obj = obj + replace( v{ i }, [ t; x{i} ], [ 0; x0{i} ] );
    end
end

opts = sdpsettings('solver','mosek');
[sol,vmono,Q,residuals,everything] = solvesos(F,-obj,opts,c);

out.sol = value(obj);
end