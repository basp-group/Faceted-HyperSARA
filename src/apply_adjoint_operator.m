function Fx = apply_adjoint_operator(y, G, varargin)
    
if isempty(varargin) % no DR
    Fx = G' * y;
else % DR
    % Sigma = varargin{1};
    Fx = G' * (varargin{1} .* y);
end 

end