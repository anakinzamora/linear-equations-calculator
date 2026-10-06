1;
% WEB_SOLVE  Runs one request from the web page through the MATLAB solver code.
%
%   octave-cli --norc --quiet --no-history web_solve.m <work folder>
%
%   Reads   <work folder>/request.json   (written by server.py)
%   Writes  <work folder>/output.txt     (the Command Window printout)
%           <work folder>/result.json    (status and the rounded unknowns)
%
%   The calculations are done by solver_core.m, which holds the same solver
%   functions as LinearSystemSolver.m and LinearSystemSolverGUI.m.

here = fileparts(mfilename('fullpath'));
source(fullfile(here, 'solver_core.m'));

args = argv();
wd = args{1};
req = jsondecode(fileread(fullfile(wd, 'request.json')));

A = double(req.A);
b = double(req.b(:));
o.dp = double(req.dp);
o.names = reshape(cellstr(req.names), 1, []);

fid = fopen(fullfile(wd, 'output.txt'), 'w');
P = @(varargin) fprintf(fid, varargin{:});
res = struct('ok', true, 'status', 'solved', 'x', {{}}, 'error', '');
try
    switch req.method
        case 'doolittle'
            x = solveDoolittle(A, b, o, P);
        case 'crout'
            x = solveCrout(A, b, o, P);
        case 'cholesky'
            x = solveCholesky(A, b, o, P);
        case 'cramer'
            x = solveCramer(A, b, o, P);
        case 'gaussjordan'
            x = solveGaussJordan(A, b, o, P);
        case {'seidel', 'jacobi'}
            it.x0 = double(req.x0(:));
            it.crit = req.crit;
            it.tol = double(req.tol);
            it.maxIt = double(req.maxIt);
            it.note = '';
            if req.reorder && ~isStrictlyDominant(A)
                [A2, b2, perm] = dominantOrder(A, b);
                if ~isempty(perm)
                    A = A2; b = b2;
                    it.note = sprintf(['The equations were rearranged into the order %s ' ...
                        'to make the system diagonally dominant.'], ...
                        strjoin(arrayfun(@(k) sprintf('(%d)', k), perm, 'UniformOutput', false), ', '));
                end
            end
            [x, res.status] = solveIterative(A, b, o, P, req.method, it);
        otherwise
            error('LSS:method', 'Unknown method.');
    end
    res.x = arrayfun(@(v) fnum(v, o.dp), x(:)', 'UniformOutput', false);
catch err
    fprintf(fid, '\n  *** CANNOT SOLVE ***\n  %s\n', err.message);
    res.ok = false;
    res.status = 'error';
    res.error = err.message;
end
fclose(fid);

fid = fopen(fullfile(wd, 'result.json'), 'w');
fprintf(fid, '%s', jsonencode(res));
fclose(fid);
