1;
% SOLVER_CORE  The solver functions of LinearSystemSolver.m / LinearSystemSolverGUI.m,
% collected in one file so GNU Octave can load them for the web version (web_solve.m).
% This file is generated from the same source as the two MATLAB files; do not edit by hand.

%% =====================================================================
%%  SOLVER CORE
%%  Every solver writes its output through P, a printf-style function
%%  handle, so the same code serves the Command Window and the GUI.
%% =====================================================================

function x = solveDoolittle(A, b, o, P)
% A.1 DOOLITTLE: [A] = [L][U] with 1s on the diagonal of [L].
    n = size(A, 1); dp = o.dp;
    L = eye(n); U = zeros(n); tol = pivotTol(A);
    hdr(P, 'A.1  DOOLITTLE METHOD  (LU Decomposition)');
    P('  [A] = [L][U], where [L] is lower triangular with 1s on its diagonal\n');
    P('  and [U] is upper triangular. Then solve [L]{Y} = {B} and [U]{X} = {Y}.\n\n');
    P('    u(i,j) = a(i,j) - SUM[k=1..i-1] l(i,k)*u(k,j)                (j >= i)\n');
    P('    l(j,i) = ( a(j,i) - SUM[k=1..i-1] l(j,k)*u(k,i) ) / u(i,i)   (j > i)\n');
    printSystem(P, A, b, o);
    printMat(P, '[A | B]', [A b], dp, n);
    sub(P, 'STEP 1: Compute the entries of [L] and [U]  (row x column)');
    for i = 1:n
        P('\n   Row %d of [U]:\n', i);
        for j = i:n
            [pS, pV, s] = pairTerms(L, U, i, j, i - 1, 'l', 'u', dp);
            U(i, j) = A(i, j) - s;
            P('     %s\n', entryLine(['u' ij(i, j)], ['a' ij(i, j)], A(i, j), pS, pV, 'plain', '', 0, U(i, j), dp));
        end
        if abs(U(i, i)) < tol
            zeroPivotError(A, ['u' ij(i, i)], 'Doolittle''s');
        end
        if i < n
            P('   Column %d of [L]:\n', i);
            for j = i + 1:n
                [pS, pV, s] = pairTerms(L, U, j, i, i - 1, 'l', 'u', dp);
                L(j, i) = (A(j, i) - s) / U(i, i);
                P('     %s\n', entryLine(['l' ij(j, i)], ['a' ij(j, i)], A(j, i), pS, pV, 'div', ['u' ij(i, i)], U(i, i), L(j, i), dp));
            end
        end
    end
    x = luFinish(P, L, U, A, b, true, false, o, '[U]');
end

function x = solveCrout(A, b, o, P)
% A.2 CROUT: [A] = [L][U] with 1s on the diagonal of [U].
    n = size(A, 1); dp = o.dp;
    L = zeros(n); U = eye(n); tol = pivotTol(A);
    hdr(P, 'A.2  CROUT''S METHOD  (LU Decomposition)');
    P('  [A] = [L][U], where [L] is lower triangular and [U] is upper\n');
    P('  triangular with 1s on its diagonal. Then solve [L]{Y} = {B} and [U]{X} = {Y}.\n\n');
    P('    l(i,j) = a(i,j) - SUM[k=1..j-1] l(i,k)*u(k,j)                (i >= j)\n');
    P('    u(j,i) = ( a(j,i) - SUM[k=1..j-1] l(j,k)*u(k,i) ) / l(j,j)   (i > j)\n');
    printSystem(P, A, b, o);
    printMat(P, '[A | B]', [A b], dp, n);
    sub(P, 'STEP 1: Compute the entries of [L] and [U]  (row x column)');
    for j = 1:n
        P('\n   Column %d of [L]:\n', j);
        for i = j:n
            [pS, pV, s] = pairTerms(L, U, i, j, j - 1, 'l', 'u', dp);
            L(i, j) = A(i, j) - s;
            P('     %s\n', entryLine(['l' ij(i, j)], ['a' ij(i, j)], A(i, j), pS, pV, 'plain', '', 0, L(i, j), dp));
        end
        if abs(L(j, j)) < tol
            zeroPivotError(A, ['l' ij(j, j)], 'Crout''s');
        end
        if j < n
            P('   Row %d of [U]:\n', j);
            for i = j + 1:n
                [pS, pV, s] = pairTerms(L, U, j, i, j - 1, 'l', 'u', dp);
                U(j, i) = (A(j, i) - s) / L(j, j);
                P('     %s\n', entryLine(['u' ij(j, i)], ['a' ij(j, i)], A(j, i), pS, pV, 'div', ['l' ij(j, j)], L(j, j), U(j, i), dp));
            end
        end
    end
    x = luFinish(P, L, U, A, b, false, true, o, '[U]');
end

function x = solveCholesky(A, b, o, P)
% A.3 CHOLESKY: [A] = [L][L]^T for a symmetric positive-definite [A].
    n = size(A, 1); dp = o.dp;
    L = zeros(n); tol = pivotTol(A);
    hdr(P, 'A.3  CHOLESKY METHOD  (LU Decomposition)');
    P('  [A] = [L][U] with [U] = [L]^T. Works only when [A] is symmetric and\n');
    P('  positive definite. Then solve [L]{Y} = {B} and [L]^T{X} = {Y}.\n\n');
    P('    l(j,j) = sqrt( a(j,j) - SUM[k=1..j-1] l(j,k)^2 )\n');
    P('    l(i,j) = ( a(i,j) - SUM[k=1..j-1] l(i,k)*l(j,k) ) / l(j,j)   (i > j)\n');
    printSystem(P, A, b, o);
    printMat(P, '[A | B]', [A b], dp, n);
    [isSym, r, c] = symmetryCheck(A);
    if ~isSym
        fail('LSS:notSymmetric', ['Cholesky''s method needs a SYMMETRIC coefficient matrix ([A] = [A]^T),\n' ...
            '  but a%s = %s while a%s = %s. Use Doolittle or Crout for this system.'], ...
            ij(r, c), fshort(A(r, c), dp), ij(c, r), fshort(A(c, r), dp));
    end
    P('\n  [A] is symmetric (a(i,j) = a(j,i)), so Cholesky''s method can be tried.\n');
    sub(P, 'STEP 1: Compute the entries of [L]  (row x column)');
    for j = 1:n
        P('\n   Column %d of [L]:\n', j);
        pS = {}; pV = {}; s = 0;
        for k = 1:j - 1
            s = s + L(j, k)^2;
            pS{end + 1} = ['l' ij(j, k) '^2']; %#ok<AGROW>
            pV{end + 1} = ['(' fnum(L(j, k), dp) ')^2']; %#ok<AGROW>
        end
        d = A(j, j) - s;
        if d <= tol
            fail('LSS:notPositiveDefinite', ['Cholesky''s method needs a POSITIVE-DEFINITE matrix, but at l%s\n' ...
                '  the value under the square root is %s (it must be greater than 0).\n' ...
                '  Use Doolittle or Crout for this system.'], ij(j, j), fnum(d, dp));
        end
        L(j, j) = sqrt(d);
        P('     %s\n', entryLine(['l' ij(j, j)], ['a' ij(j, j)], A(j, j), pS, pV, 'sqrt', '', 0, L(j, j), dp));
        for i = j + 1:n
            pS = {}; pV = {}; s = 0;
            for k = 1:j - 1
                s = s + L(i, k) * L(j, k);
                pS{end + 1} = ['l' ij(i, k) '*l' ij(j, k)]; %#ok<AGROW>
                pV{end + 1} = ['(' fnum(L(i, k), dp) ')(' fnum(L(j, k), dp) ')']; %#ok<AGROW>
            end
            L(i, j) = (A(i, j) - s) / L(j, j);
            P('     %s\n', entryLine(['l' ij(i, j)], ['a' ij(i, j)], A(i, j), pS, pV, 'div', ['l' ij(j, j)], L(j, j), L(i, j), dp));
        end
    end
    x = luFinish(P, L, L', A, b, false, false, o, '[U] = [L]^T');
end

function x = luFinish(P, L, U, A, b, unitL, unitU, o, uLabel)
% Shared ending of the three LU methods: show [L], [U], {Y} and {X}.
    dp = o.dp; n = numel(b);
    sub(P, 'STEP 2: The [L] and [U] matrices');
    printMat(P, '[L]', L, dp, 0);
    printMat(P, uLabel, U, dp, 0);
    P('\n  Check: max |[L][U] - [A]| = %.2e\n', max(max(abs(L * U - A))));
    sub(P, 'STEP 3: Forward substitution  [L]{Y} = {B}');
    y = forwardSub(P, L, b, unitL, dp);
    yNames = arrayfun(@(k) sprintf('y%d', k), 1:n, 'UniformOutput', false);
    printVec(P, '{Y}  (individual variable matrix)', y, yNames, dp);
    sub(P, ['STEP 4: Back substitution  ' uLabel(1:3) '{X} = {Y}']);
    x = backSub(P, U, y, unitU, o);
    printVec(P, '{X}  (the unknowns)', x, o.names, dp);
    printAnswer(P, x, A, b, o);
end

function y = forwardSub(P, L, b, unitDiag, dp)
    n = numel(b); y = zeros(n, 1);
    for i = 1:n
        pS = {}; pV = {}; s = 0;
        for k = 1:i - 1
            s = s + L(i, k) * y(k);
            pS{end + 1} = sprintf('l%s*y%d', ij(i, k), k); %#ok<AGROW>
            pV{end + 1} = sprintf('(%s)(%s)', fnum(L(i, k), dp), fnum(y(k), dp)); %#ok<AGROW>
        end
        if unitDiag
            y(i) = b(i) - s;
            line = entryLine(sprintf('y%d', i), sprintf('b%d', i), b(i), pS, pV, 'plain', '', 0, y(i), dp);
        else
            y(i) = (b(i) - s) / L(i, i);
            line = entryLine(sprintf('y%d', i), sprintf('b%d', i), b(i), pS, pV, 'div', ['l' ij(i, i)], L(i, i), y(i), dp);
        end
        P('     %s\n', line);
    end
end

function x = backSub(P, U, y, unitDiag, o)
    n = numel(y); x = zeros(n, 1); dp = o.dp; nm = o.names;
    for i = n:-1:1
        pS = {}; pV = {}; s = 0;
        for k = i + 1:n
            s = s + U(i, k) * x(k);
            pS{end + 1} = sprintf('u%s*%s', ij(i, k), nm{k}); %#ok<AGROW>
            pV{end + 1} = sprintf('(%s)(%s)', fnum(U(i, k), dp), fnum(x(k), dp)); %#ok<AGROW>
        end
        if unitDiag
            x(i) = y(i) - s;
            line = entryLine(nm{i}, sprintf('y%d', i), y(i), pS, pV, 'plain', '', 0, x(i), dp);
        else
            x(i) = (y(i) - s) / U(i, i);
            line = entryLine(nm{i}, sprintf('y%d', i), y(i), pS, pV, 'div', ['u' ij(i, i)], U(i, i), x(i), dp);
        end
        P('     %s\n', line);
    end
end

function x = solveCramer(A, b, o, P)
% B. CRAMER'S RULE: x(i) = det[A_i] / det[A].
    n = size(A, 1); dp = o.dp; nm = o.names;
    hdr(P, 'B.  CRAMER''S RULE');
    P('  x(i) = det[A_i] / det[A], where [A_i] is [A] with column i\n');
    P('  replaced by the constants {B}.\n');
    printSystem(P, A, b, o);
    sub(P, 'STEP 1: The main determinant, det[A]');
    printMat(P, '[A]', A, dp, 0);
    printMat(P, '{B}', b, dp, 0);
    D = cleanDet(det(A), A);
    P('\n     det[A] = %s\n', fnum(D, dp));
    if D == 0 || rcond(A) < eps
        fail('LSS:singular', ['det[A] = 0, so the system has no unique solution\n' ...
            '  (the equations are dependent or inconsistent). Cramer''s rule cannot be applied.']);
    end
    sub(P, 'STEP 2: The matrices and determinants of the unknowns');
    Ds = zeros(n, 1);
    for i = 1:n
        Ai = A; Ai(:, i) = b;
        Ds(i) = cleanDet(det(Ai), Ai);
        P('\n  (%d) Unknown %s: replace column %d of [A] with {B}\n', i, nm{i}, i);
        printMat(P, sprintf('[A_%s]', nm{i}), Ai, dp, 0);
        P('     det[A_%s] = %s\n', nm{i}, fnum(Ds(i), dp));
    end
    sub(P, 'STEP 3: The values of the unknowns');
    x = Ds / D;
    for i = 1:n
        P('     %s = det[A_%s] / det[A] = %s / %s = %s\n', nm{i}, nm{i}, ...
            fnum(Ds(i), dp), fnum(D, dp), fnum(x(i), dp));
    end
    printAnswer(P, x, A, b, o);
end

function x = solveGaussJordan(A, b, o, P)
% C. GAUSS-JORDAN: reduce [A | B] to [I | X] with elementary row operations.
    n = size(A, 1); dp = o.dp; nm = o.names;
    M = [A b]; tol = pivotTol(A); stepNo = 0;
    hdr(P, 'C.  GAUSS-JORDAN ELIMINATION');
    P('  Goal: use row operations to turn the augmented matrix [A | B] into\n');
    P('  [I | X]. Each column is processed in turn: make the pivot 1, then\n');
    P('  make every other entry in that column 0.\n');
    printSystem(P, A, b, o);
    printMat(P, 'Initial augmented matrix [A | B]', M, dp, n);
    for k = 1:n
        sub(P, sprintf('COLUMN %d  (pivot = row %d, column %d)', k, k, k));
        if abs(M(k, k)) < tol
            [~, rel] = max(abs(M(k:n, k)));
            p = k - 1 + rel;
            if abs(M(p, k)) < tol
                fail('LSS:singular', ['Column %d has no non-zero pivot, so the matrix is singular and\n' ...
                    '  the system has no unique solution.'], k);
            end
            M([k p], :) = M([p k], :);
            stepNo = stepNo + 1;
            P('\n  Step %d:  R%d <-> R%d   (the pivot is 0, so swap rows)\n', stepNo, k, p);
            printMat(P, 'Result', M, dp, n);
        end
        piv = M(k, k);
        if abs(piv - 1) > 1e-12
            M(k, :) = M(k, :) / piv;
            M(k, k) = 1;
            M(k, abs(M(k, :)) < 1e-13) = 0;
            stepNo = stepNo + 1;
            P('\n  Step %d:  R%d <- R%d / (%s)   (make the pivot equal to 1)\n', stepNo, k, k, fnum(piv, dp));
            printMat(P, 'Result', M, dp, n);
        else
            P('\n  The pivot is already 1.\n');
        end
        for i = [1:k - 1, k + 1:n]
            f = M(i, k);
            if f == 0, continue; end
            M(i, :) = M(i, :) - f * M(k, :);
            M(i, k) = 0;
            M(i, abs(M(i, :)) < 1e-13 * max(1, max(abs(M(i, :))))) = 0;
            stepNo = stepNo + 1;
            if f > 0
                op = sprintf('R%d <- R%d - (%s)R%d', i, i, fnum(f, dp), k);
            else
                op = sprintf('R%d <- R%d + (%s)R%d', i, i, fnum(-f, dp), k);
            end
            P('\n  Step %d:  %s   (make a%s = 0)\n', stepNo, op, ij(i, k));
            printMat(P, 'Result', M, dp, n);
        end
    end
    sub(P, 'FINAL MATRIX  (before the final substitution)');
    printMat(P, '[I | X]', M, dp, n);
    sub(P, 'FINAL SUBSTITUTION  (read each row of [I | X])');
    x = M(:, n + 1);
    for i = 1:n
        P('     Row %d:  (1)%s = %s   ->   %s = %s\n', i, nm{i}, fnum(x(i), dp), nm{i}, fnum(x(i), dp));
    end
    printAnswer(P, x, A, b, o);
end

function [x, status] = solveIterative(A, b, o, P, method, it)
% D. GAUSS-SEIDEL ('seidel') or GAUSS-JACOBI ('jacobi').
%   it.x0 (initial guess), it.crit ('abs' | 'rel' | 'repeat'),
%   it.tol, it.maxIt, it.note (text shown when the rows were rearranged)
    n = size(A, 1); dp = o.dp; nm = o.names;
    isSeidel = strcmp(method, 'seidel');
    if isSeidel
        hdr(P, 'D.1  GAUSS-SEIDEL METHOD  (Iterative)');
        P('  Each equation is solved for its diagonal unknown. Gauss-Seidel uses\n');
        P('  every NEW value as soon as it is computed in the same iteration.\n');
    else
        hdr(P, 'D.2  GAUSS-JACOBI METHOD  (Iterative)');
        P('  Each equation is solved for its diagonal unknown. Gauss-Jacobi uses\n');
        P('  only the values from the PREVIOUS iteration to compute new ones.\n');
    end
    if isfield(it, 'note') && ~isempty(it.note)
        P('\n  NOTE: %s\n', it.note);
    end
    printSystem(P, A, b, o);
    if rcond(A) < eps
        fail('LSS:singular', ['det[A] = 0, so the system has no unique solution\n' ...
            '  (the equations are dependent or inconsistent).']);
    end
    for i = 1:n
        if A(i, i) == 0
            fail('LSS:zeroDiagonal', ['a%s = 0. Iterative methods divide by the diagonal entries, so every\n' ...
                '  a(i,i) must be non-zero. Rearrange the equations and try again.'], ij(i, i));
        end
    end

    sub(P, 'STEP 1: Check for diagonal dominance');
    strict = true;
    for i = 1:n
        d = abs(A(i, i)); offSum = sum(abs(A(i, :))) - d;
        if d > offSum, verdict = '>   OK'; else, verdict = '<=  not dominant'; strict = false; end
        others = '';
        for j = [1:i - 1, i + 1:n]
            if isempty(others), others = sprintf('|a%s|', ij(i, j));
            else, others = sprintf('%s + |a%s|', others, ij(i, j)); end
        end
        P('     Row %d:  |a%s| = %s   vs   %s = %s   ->  %s\n', i, ij(i, i), fnum(d, dp), ...
            others, fnum(offSum, dp), verdict);
    end
    if strict
        P('\n  The system is strictly diagonally dominant, so the method will converge.\n');
    else
        P('\n  WARNING: The system is NOT strictly diagonally dominant, so convergence\n');
        P('  is not guaranteed. The iterations may diverge.\n');
    end

    sub(P, 'STEP 2: Iteration formulas');
    for i = 1:n
        f = sprintf('%s(k+1) = ( %s', nm{i}, fnum(b(i), dp));
        for j = [1:i - 1, i + 1:n]
            if isSeidel && j < i, kk = 'k+1'; else, kk = 'k'; end
            c = -A(i, j);
            if c == 0, continue; end
            if c < 0, sg = '-'; else, sg = '+'; end
            f = sprintf('%s %s %s*%s(%s)', f, sg, fnum(abs(c), dp), nm{j}, kk);
        end
        P('     %s ) / %s\n', f, fnum(A(i, i), dp));
    end

    x0 = it.x0(:);
    switch it.crit
        case 'abs'
            critText = sprintf('stop when every |x(k) - x(k-1)| <= %s', fshortTol(it.tol));
            errLabel = 'Ea';
        case 'rel'
            critText = sprintf('stop when every |(x(k) - x(k-1)) / x(k)| x 100%% <= %s%%', fshortTol(it.tol));
            errLabel = 'Ea%';
        otherwise
            critText = sprintf('stop when the values repeat to %d decimal places', dp);
            errLabel = 'Ea';
    end
    P('\n  Initial guess:      %s\n', vecInline(x0, nm, dp));
    P('  Stopping criterion: %s\n', critText);
    P('  Maximum iterations: %d\n', it.maxIt);

    X = zeros(it.maxIt + 1, n); E = nan(it.maxIt + 1, n);
    X(1, :) = x0';
    xo = x0; status = 'maxit'; K = it.maxIt;
    limit = 1e6 * max([1; abs(b(:)); abs(x0(:))]);
    for k = 1:it.maxIt
        xn = xo;
        for i = 1:n
            if isSeidel, src = xn; else, src = xo; end
            s = b(i) - A(i, :) * src + A(i, i) * src(i);
            xn(i) = s / A(i, i);
        end
        switch it.crit
            case 'rel'
                e = abs((xn - xo) ./ xn) * 100;
                e(xn == 0 & xo == 0) = 0;
                e(xn == 0 & xo ~= 0) = Inf;
                done = all(e <= it.tol);
            case 'repeat'
                e = abs(xn - xo);
                done = all(rhe(xn, dp) == rhe(xo, dp));
            otherwise
                e = abs(xn - xo);
                done = all(e <= it.tol);
        end
        X(k + 1, :) = xn'; E(k + 1, :) = e';
        if any(~isfinite(xn)) || max(abs(xn)) > limit
            status = 'diverged'; K = k; xo = xn; break;
        end
        xo = xn;
        if done
            status = 'converged'; K = k; break;
        end
    end
    x = xo;

    sub(P, sprintf('STEP 3: Iterations  (%d shown)', K));
    printIterTable(P, X(1:K + 1, :), E(1:K + 1, :), nm, dp, errLabel);
    switch status
        case 'converged'
            P('\n  CONVERGED after %d iteration(s): the stopping criterion is satisfied.\n', K);
            printAnswer(P, x, A, b, o);
        case 'diverged'
            P('\n  DIVERGED: the values grow without bound, so this arrangement of the\n');
            P('  equations does not converge. Rearrange the equations so the system is\n');
            P('  diagonally dominant, or use a direct method (A, B or C).\n');
        otherwise
            P('\n  STOPPED after the maximum of %d iterations without meeting the\n', it.maxIt);
            P('  stopping criterion. The last values are shown below.\n');
            printAnswer(P, x, A, b, o);
    end
end

function [A2, b2, perm] = dominantOrder(A, b)
% Row order that makes [A] strictly diagonally dominant ([] if none exists).
    n = size(A, 1); perm = zeros(1, n); A2 = A; b2 = b;
    for i = 1:n
        [m, c] = max(abs(A(i, :)));
        if m <= sum(abs(A(i, :))) - m || perm(c) ~= 0
            perm = []; return;
        end
        perm(c) = i;
    end
    A2 = A(perm, :); b2 = b(perm);
end

function tf = isStrictlyDominant(A)
    d = abs(diag(A));
    tf = all(d > sum(abs(A), 2) - d);
end

%% ---------------------------------------------------------------------
%%  Output helpers
%% ---------------------------------------------------------------------

function v = rhe(x, dp)
% Round to dp decimal places, ties to even (banker's rounding).
    s = 10^dp;
    y = x .* s;
    v = round(y);
    fr = abs(y - fix(y));
    tie = abs(fr - 0.5) <= max(1e-9, abs(y) .* 1e-12);
    v(tie) = 2 .* round(y(tie) ./ 2);
    v = v ./ s;
    v(v == 0) = 0;
end

function s = fnum(x, dp)
% Number rounded to dp decimals (ties to even), as text.
    if isnan(x), s = 'NaN'; return; end
    if isinf(x)
        if x > 0, s = 'Inf'; else, s = '-Inf'; end
        return;
    end
    v = rhe(x, dp);
    if abs(v) >= 1e10
        s = sprintf(['%.' num2str(dp) 'e'], x);
    else
        s = sprintf(['%.' num2str(dp) 'f'], v);
    end
end

function s = fshort(x, dp)
% Like fnum, without trailing zeros (used inside equations).
    s = fnum(x, dp);
    if any(s == '.') && ~any(s == 'e')
        s = regexprep(s, '0+$', '');
        s = regexprep(s, '\.$', '');
    end
end

function s = fshortTol(t)
    s = regexprep(sprintf('%.10f', t), '0+$', '');
    s = regexprep(s, '\.$', '');
end

function s = ij(i, j)
% Subscript text: 12 for (1,2); (10,11) once an index reaches 10.
    if i < 10 && j < 10
        s = sprintf('%d%d', i, j);
    else
        s = sprintf('(%d,%d)', i, j);
    end
end

function zeroPivotError(A, entry, method)
    if rcond(A) < eps
        fail('LSS:singular', ['Zero pivot: %s = 0 because [A] is singular (det[A] = 0).\n' ...
            '  The system has no unique solution.'], entry);
    end
    fail('LSS:zeroPivot', ['Zero pivot: %s = 0, so %s method (no row exchanges) cannot continue.\n' ...
        '  Rearrange the equations (swap rows) and try again, or use Gauss-Jordan.'], entry, method);
end

function fail(id, fmt, varargin)
% Throw an error whose message is fully formatted (works the same in every MATLAB release).
    error(id, '%s', sprintf(fmt, varargin{:}));
end

function tol = pivotTol(A)
    tol = 1e-12 * max(1, max(abs(A(:))));
end

function D = cleanDet(D, M)
% Integer matrices have integer determinants: remove floating-point noise.
    if all(M(:) == round(M(:))) && abs(D - round(D)) < 1e-6 * max(1, abs(D))
        D = round(D);
    end
end

function [tf, r, c] = symmetryCheck(A)
    tf = true; r = 0; c = 0;
    n = size(A, 1); tol = 1e-10 * max(1, max(abs(A(:))));
    for i = 1:n
        for j = i + 1:n
            if abs(A(i, j) - A(j, i)) > tol
                tf = false; r = i; c = j; return;
            end
        end
    end
end

function [pS, pV, s] = pairTerms(L, U, r, c, m, ln, un, dp)
% Terms l(r,k)*u(k,c) for k = 1..m, as symbols, substituted values and their sum.
    pS = {}; pV = {}; s = 0;
    for k = 1:m
        s = s + L(r, k) * U(k, c);
        pS{end + 1} = [ln ij(r, k) '*' un ij(k, c)]; %#ok<AGROW>
        pV{end + 1} = ['(' fnum(L(r, k), dp) ')(' fnum(U(k, c), dp) ')']; %#ok<AGROW>
    end
end

function line = entryLine(name, aSym, aVal, pS, pV, mode, dSym, dVal, res, dp)
% Builds  name = formula = substituted values = result
%   mode 'plain': a - SUM,  'div': (a - SUM) / d,  'sqrt': sqrt(a - SUM)
    sym = aSym; num = fnum(aVal, dp);
    for k = 1:numel(pS)
        sym = [sym ' - ' pS{k}]; %#ok<AGROW>
        num = [num ' - ' pV{k}]; %#ok<AGROW>
    end
    hasTerms = ~isempty(pS);
    switch mode
        case 'div'
            if hasTerms
                sym = ['(' sym ') / ' dSym];
                num = ['(' num ') / ' fnum(dVal, dp)];
            else
                sym = [sym ' / ' dSym];
                num = [num ' / ' fnum(dVal, dp)];
            end
        case 'sqrt'
            sym = ['sqrt(' sym ')'];
            num = ['sqrt(' num ')'];
    end
    if strcmp(mode, 'plain') && ~hasTerms
        line = sprintf('%s = %s = %s', name, sym, fnum(res, dp));
    else
        line = sprintf('%s = %s = %s = %s', name, sym, num, fnum(res, dp));
    end
end

function hdr(P, t)
    bar = repmat('=', 1, 70);
    P('\n%s\n  %s\n%s\n', bar, t, bar);
end

function sub(P, t)
    P('\n%s\n  %s\n', repmat('-', 1, 70), t);
end

function S = fmtCells(M, dp)
% Every entry of M as text (same result as fnum), formatted in one sprintf call.
    [r, c] = size(M);
    V = rhe(M, dp);
    if any(abs(V(:)) >= 1e10 & isfinite(V(:)))
        S = cell(r, c);
        for i = 1:r
            for j = 1:c
                S{i, j} = fnum(M(i, j), dp);
            end
        end
        return;
    end
    t = strsplit(sprintf(['%.' num2str(dp) 'f\n'], V.'), sprintf('\n'));
    S = reshape(t(1:r * c), c, r).';
end

function lines = matLines(M, dp, nAug)
    S = fmtCells(M, dp);
    [r, c] = size(S);
    w = max(cellfun('length', S(:)));
    cellFmt = ['  %' num2str(w) 's'];
    lines = cell(r, 1);
    for i = 1:r
        if nAug > 0 && nAug < c
            lines{i} = ['[' sprintf(cellFmt, S{i, 1:nAug}) '  |' sprintf(cellFmt, S{i, nAug + 1:c}) '  ]'];
        else
            lines{i} = ['[' sprintf(cellFmt, S{i, :}) '  ]'];
        end
    end
end

function printMat(P, label, M, dp, nAug)
    P('\n  %s =\n', label);
    L = matLines(M, dp, nAug);
    for i = 1:numel(L)
        P('     %s\n', L{i});
    end
end

function printVec(P, label, v, rowNames, dp)
    P('\n  %s =\n', label);
    L = matLines(v(:), dp, 0);
    for i = 1:numel(L)
        P('     %s   %s\n', L{i}, rowNames{i});
    end
end

function s = vecInline(v, names, dp)
    parts = cell(1, numel(v));
    for i = 1:numel(v)
        parts{i} = sprintf('%s = %s', names{i}, fnum(v(i), dp));
    end
    s = strjoin(parts, ',  ');
end

function s = eqStr(a, bi, names, dp)
% One equation as text, e.g.  4x1 - x2 + x3 = 12
    s = '';
    for j = 1:numel(a)
        c = a(j);
        if c == 0, continue; end
        if isempty(s)
            if c < 0, sg = '-'; else, sg = ''; end
        else
            if c < 0, sg = ' - '; else, sg = ' + '; end
        end
        if abs(c) == 1, cs = ''; else, cs = fshort(abs(c), dp); end
        s = [s sg cs names{j}]; %#ok<AGROW>
    end
    if isempty(s), s = '0'; end
    s = [s ' = ' fshort(bi, dp)];
end

function printSystem(P, A, b, o)
    n = size(A, 1);
    P('\n  System of %d linear equations:\n', n);
    for i = 1:n
        P('     (%d)  %s\n', i, eqStr(A(i, :), b(i), o.names, o.dp));
    end
end

function printAnswer(P, x, A, b, o)
    sub(P, 'ANSWER');
    for i = 1:numel(x)
        P('     %s = %s\n', o.names{i}, fnum(x(i), o.dp));
    end
    P('\n  (Values are rounded to %d decimal places, ties to even.)\n', o.dp);
    P('  Check: max |[A]{X} - {B}| = %.2e\n', max(abs(A * x(:) - b(:))));
end

function printIterTable(P, X, E, names, dp, errLabel)
    [m, n] = size(X);
    heads = cell(1, 1 + 2 * n);
    heads{1} = 'k';
    for j = 1:n
        heads{1 + j} = names{j};
        heads{1 + n + j} = sprintf('%s(%s)', errLabel, names{j});
    end
    C = [arrayfun(@(k) sprintf('%d', k), (0:m - 1)', 'UniformOutput', false), ...
        fmtCells(X, dp), fmtCells(E, dp)];
    C(1, n + 2:end) = {'-'};
    w = zeros(1, 1 + 2 * n);
    for j = 1:numel(w)
        w(j) = max([length(heads{j}), cellfun(@length, C(:, j))']);
    end
    line = ' ';
    for j = 1:numel(w)
        line = [line '  ' sprintf(['%' num2str(w(j)) 's'], heads{j})]; %#ok<AGROW>
        if j == 1 + n, line = [line '  |']; end %#ok<AGROW>
    end
    P('\n  %s\n', line);
    P('  %s\n', repmat('-', 1, length(line)));
    for r = 1:m
        line = ' ';
        for j = 1:numel(w)
            line = [line '  ' sprintf(['%' num2str(w(j)) 's'], C{r, j})]; %#ok<AGROW>
            if j == 1 + n, line = [line '  |']; end %#ok<AGROW>
        end
        P('  %s\n', line);
    end
end

function nm = defaultNames(n)
    nm = arrayfun(@(k) sprintf('x%d', k), 1:n, 'UniformOutput', false);
end

function [A, b] = exampleSystem(n)
% Symmetric, positive-definite and diagonally dominant, so every method works.
    if n == 3
        A = [4 -1 1; -1 4 -2; 1 -2 4];
        b = [12; -1; 5];
    else
        A = 4 * eye(n) - diag(ones(n - 1, 1), 1) - diag(ones(n - 1, 1), -1);
        b = A * (1:n)';
    end
end
