function y = apply_direct_operator(Fx, G, Sigma)

y = Sigma .* (G * Fx);

end