function [u,v,w,na] = generate_uv_cov_antennas(Pos_ant,x0,h,lat,dec, nstep)


% creation positions u-v des antennes
na = max(size(Pos_ant)) ; % nombre d antennes

% antennes originales
b = zeros(na,2) ;
b(:,1) = sqrt( (Pos_ant(:,1)-x0(1)).^2 + (Pos_ant(:,2)-x0(2)).^2 ) ;
b(:,2) = atan2( Pos_ant(:,1)-x0(1) , Pos_ant(:,2)-x0(2) ) ;

% selection dans le plan u-v
u=cell(na,1);
v=cell(na,1);
w=cell(na,1);
for alpha = 1:na
X = uvTrack(h, b(alpha, 1), b(alpha, 2), 0., lat, dec, nstep);
u{alpha}=X(:,1);
v{alpha}=X(:,2);
w{alpha}=X(:,3);
end




% % % % selection dans le plan u-v
% % % u=[];
% % % v=[];
% % % w=[];
% % % for alpha = 1:na
% % % X = uvTrack(h, b(alpha, 1), b(alpha, 2), 0., lat, dec, nstep);
% % % u=[u ; X(:,1)];
% % % v=[v ; X(:,2)];
% % % w=[w ; X(:,3)];
% % % end

U = cell2mat(u) ;
V = cell2mat(v) ;

% % % figure
% % % plot(U, V, '.')
% % % hold on
% % % plot(-U, -V, '.r')
% % % xlabel('cont u-v cov')

end