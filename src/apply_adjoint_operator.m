function Fx = apply_adjoint_operator(y, G, Sigma)
    
% if isempty(varargin) % no DR
%     Fx = G' * y;
% else % DR
%     % Sigma = varargin{1};
Fx = G' * (Sigma .* y);
% end 

end