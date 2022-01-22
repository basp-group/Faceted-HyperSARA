function [Anorm, squared_operator_norm, rel_var, squared_precond_operator_norm, rel_var_precond] = util_operator_norm(G, W, A, At, aW, Ny, Nx, tol, max_iter,flag_dr,Sigma)
% [extended_summary]
% 
% Parameters
% ----------
% G : [type]
%     [description]
% W : [type]
%     [description]
% A : [type]
%     [description]
% At : [type]
%     [description]
% aW : [type]
%     [description]
% Ny : [type]
%     [description]
% Nx : [type]
%     [description]
% tol : [type]
%     [description]
% max_iter : [type]
%     [description]
% flag_dr : [type]
%     [description]
% Sigma : [type]
%     [description]
% 
% Returns
% -------
% [type]
%     [description]

if nargin<10 
    flag_dr=0; Sigma =[];
end
nchans = numel(G);
squared_precond_operator_norm = zeros(nchans, 1);
rel_var_precond = zeros(nchans, 1);

for l = 1:nchans
    F = afclean(@(x) HS_forward_operator_precond_G(x, G(l), W(l), A, aW(l),flag_dr,Sigma));
    Ft = afclean(@(y) HS_adjoint_operator_precond_G(y, G(l), W(l), At, aW(l), Ny, Nx,flag_dr,Sigma));
    [squared_precond_operator_norm(l), rel_var_precond(l)] = op_norm(F, Ft, [Ny, Nx], tol, max_iter, 0);
end

% ! operator is block diagonal: it norm is thus given by the max singular value
% ! make sure the pdfb convergence criterion is strictly satisfied
% ! by taking a smaller step-size (take precision of the estimation
% ! of the operator norm into account)
Anorm = max(squared_precond_operator_norm .* (1 + rel_var_precond));

squared_operator_norm = zeros(nchans, 1);
rel_var = zeros(nchans, 1);
for l = 1:nchans
    F = afclean(@(x) HS_forward_operator_G(x, G(l), W(l), A,flag_dr,Sigma));
    Ft = afclean(@(y) HS_adjoint_operator_G(y, G(l), W(l), At, Ny, Nx,flag_dr,Sigma));
    [squared_operator_norm(l), ~] = op_norm(F, Ft, [Ny, Nx], tol, max_iter, 0);
end

end
