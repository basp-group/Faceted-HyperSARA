function x_overlap = comm2d_update_ghost_cells_wideband(x_overlap, overlap_q, overlap_g_south_east, overlap_g_south, overlap_g_east, Qy, Qx, Qc)
% comm2d_update_ghost_cells: update the ghost cells (overlap regions) in a 
% 2D image tessellation (2D communication grid).
%-------------------------------------------------------------------------%
%%
% Input:
%
% > x_overlap             spatial facet considered (with ghost cells)
% > overlap_q             size of the left and top facet extensions 
%                         (ghost cells) [2 ,1] 
% > overlap_g_south_east  size of the ghost cells for the south-east
%                         neighbour [2, 1]
% > overlap_g_south       size of the ghost cells for the southern
%                         neighbour [2, 1]
% > overlap_g_east        size of the ghost cells for the eastern
%                         neighbour [2, 1]
% > Qy                    number of spatial facets along the y axis [1]
% > Qx                    number of spatial facets along the x axis [1]
% > Qc                    number of spectral facets along the z axis [1]
%
% Output:
%
% < x_overlap             updated spatial facet
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% [02/04/2019] final debug, code ok
% [03/04/2019] add basic support for wideband, see if further modifications 
% are needed later on (for the spectral splitting)
%-------------------------------------------------------------------------%
%%

[i, q] = ind2sub([Qc, Qy*Qx], labindex);
[qy, qx] = ind2sub([Qy, Qx], q);
get_index = @(i, q) (q-1)*Qc + i;

% destination
dest_vert = []; % workers sending information to the current worker
dest_horz = []; % workers waiting for informations from the current worker
dest_diag = [];

% reception
src_vert = []; % workers sending information to the current worker
src_horz = []; % workers waiting for informations from the current worker
src_diag = []; % up to 3 neighbours for this communication, to be defined later on

% data to send
data_vert = [];
data_horz = [];
data_diag = [];

% define communications (to S, E, SE)
if qy < Qy
    % S communication (qy+1, qx) -> q = (qx-1)*Qy + qy+1
    dest_vert = get_index(i, (qx-1)*Qy + qy+1);
    if qx > 1
        offset_x = overlap_g_south(2);
    else
        offset_x = 0;
    end
    data_vert = x_overlap(end-overlap_g_south(1)+1:end, 1+offset_x:end, :);
    %data_vert = labindex*ones(size(data_vert)); % just for testing purposes
    
    if qx < Qx
        % communicate corner to SE (qy+1, qx+1) -> q = qx*Qy + qy+1
        dest_diag = get_index(i, qx*Qy + qy+1);
        data_diag = x_overlap(end-overlap_g_south_east(1)+1:end, end-overlap_g_south_east(2)+1:end, :);
        %data_diag = labindex*ones(size(data_diag)); % just for testing purposes
    end
end

if qx < Qx
    % E communication (qy, qx+1) -> q = qx*Qy + qy
    dest_horz = get_index(i, qx*Qy + qy);
    if qy > 1
        offset_y = overlap_g_east(1);
    else
        offset_y = 0;
    end
    data_horz = x_overlap(1+offset_y:end, end-overlap_g_east(2)+1:end, :);
    %data_horz = labindex*ones(size(data_horz)); % just for testing purposes
end

% define receptions (from N, W, NW)
if qy > 1
    % N reception (qy-1, qx) -> q = (qx-1)*Qy + qy-1
    src_vert = get_index(i, (qx-1)*Qy + qy-1);
    if qx > 1
        % NW reception (qy-1, qx-1) -> q = (qx-2)*Qy + qy-1
        src_diag = get_index(i, (qx-2)*Qy + qy-1);
    end
end

if qx > 1
    % W reception (qy, qx-1) -> q = (qx-2)*Qy + qy
    src_horz = get_index(i, (qx-2)*Qy + qy);
end

% is there a way to do everything at once? (i.e., reduce the
% synchronization phasea induced by the three calls to labSendReceive?)
% see if this can be as efficient as the gop (to be confirmed)
% vertical communications
rcv_vert_data = labSendReceive(dest_vert, src_vert, data_vert);
% horizontal communications
rcv_horz_data = labSendReceive(dest_horz, src_horz, data_horz);
% diagonal communications
rcv_diag_data = labSendReceive(dest_diag, src_diag, data_diag);


% update portions of the overlapping facet with the received data
if ~isempty(rcv_vert_data) % from N
    % offset if qx > 1
    if qx > 1
        offset = overlap_q(2);
    else
        offset = 0;
    end
    % x_overlap(1:size(rcv_vert_data, 1), 1+size(x_overlap,2)-size(rcv_vert_data,2):end) = rcv_vert_data;
    x_overlap(1:size(rcv_vert_data, 1), 1+offset:end, :) = rcv_vert_data;
end

% from W
if ~isempty(rcv_horz_data)
    % offset qy > 1
    if qy > 1
        offset = overlap_q(1);
    else
        offset = 0;
    end
    %x_overlap(1+offset:end, 1:overlap_q(2)) = rcv_horz_data;
    x_overlap(1+offset:end, 1:size(rcv_horz_data, 2), :) = rcv_horz_data;
end

% from NW
if ~isempty(rcv_diag_data)
    x_overlap(1:size(rcv_diag_data, 1), 1:size(rcv_diag_data, 2), :) = rcv_diag_data;
end

end
