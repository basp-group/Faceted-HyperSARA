function [u_ab, v_ab, na, antennas] = generate_uv_coverage(T, hrs, dl, cov_type, varargin)

% ----------------------------------------------------- %
switch cov_type                              
    % VLA coverage ---------------------------------------- %
    case 'vlaa'                                             
        fileID = fopen('vlaa.itrf-copie.txt') ;                 
        save_cov_file = ['gen_uv_tracks/ant_vla_pos.mat'];      
    % ASKA coverage --------------------------------------- %
    case 'askap'                                            
        fileID = fopen('askap.itrf-copie.txt') ;                
        save_cov_file = ['gen_uv_tracks/ant_askap_pos.mat'];    
    % MeerKAT coverage ------------------------------------ %
    case 'meerkat'                                         
        fileID = fopen('MeerKAT.enu.txt') ;                     
        save_cov_file = ['gen_uv_tracks/ant_meerkat_pos.mat'];  
    % random cont. antenna positions                        %
    case 'random'        
        na = varargin{1};
        uv(1,:) = rand(1,2) ;                                   
        for alpha = 2:na                             
            uv_ = rand(1,2) ;                                       
            while ismember(uv_,uv,'rows')                           
                uv_ = rand(1,2) ;                                      
            end                                                     
            uv(alpha,:) = uv_ ;                                     
        end                                                    
        antenna_position = 1e06*[uv, zeros(na,1)] ;           
        save_cov_file = ['gen_uv_tracks/rand_pos.mat'];         
end                                                     
% ----------------------------------------------------- %
if strcmp(cov_type, 'rand_ant') == 0
    C = textscan(fileID,'%f %f %f %s %s %s');
    antenna_position = cell2mat({C{1} C{2} C{3}}) ;
    na = max(size(antenna_position)) ;
    fclose(fileID);
end

% M = T*na*(na-1)/2 ;
% ant_pair = na*(na-1)/2 ;

% Physical parameters
h = linspace(-hrs, hrs, T)*pi/12;               % hour angle range of +/-  hours  5
dec = (pi/180)*(40);                        % Cas A is 56.4
lat = (pi/180)*(38. + 59./60. + 48./3600.); % College Park
% reference position of x0 in the x-y plane
x0 = [mean(antenna_position(:,1)), mean(antenna_position(:,2))];
[u, v, ~] = generate_uv_cov_antennas(antenna_position, x0, h, lat, dec, T); % [u, v, w], [T, na]

M = na*(na-1)/2;
u_ab = zeros(T, M);
v_ab = zeros(T, M);
q = 1;
for a = 1:na-1
    for b = a+1:na
        u_ab(:,q) = u(:,a) - u(:,b); % u_alpha - u_beta
        v_ab(:,q) = v(:,a) - v(:,b); % v_alpha - v_beta
        q = q+1;
    end
end
            
% Antenna positions (u_ab, v_ab)
bmax = max(sqrt(u_ab(:).^2 + v_ab(:).^2));

% setting imaging params
% dl = pixel_size; % ask Arwa
v_ab = v_ab.'*pi/(bmax*dl); % dl = 1.1
u_ab = u_ab.'*pi/(bmax*dl);

% generate corresponding antenna indices (couple)
Ant_T = zeros(M, 2);
for a = 1:na
    Ant_T((a-1)*na-(a-1)*a/2+1:a*na-a*(a+1)/2, 1) = a*ones((na-a), 1); % antenna alpha
    Ant_T((a-1)*na-(a-1)*a/2+1:a*na-a*(a+1)/2, 2) = (a+1:na).';        % antenna beta
end
antennas = repmat(Ant_T, [T, 1]); % see if couples (beta, alpha are necessary)

end