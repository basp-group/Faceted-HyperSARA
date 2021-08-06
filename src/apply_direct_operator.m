function y = apply_direct_operator(Fx, G, Sigma)

% if isempty(varargin) % no DR
%     y = G * Fx;
% else % DR
%     % Sigma = varargin{1};
y = Sigma .* (G * Fx);
% end

end