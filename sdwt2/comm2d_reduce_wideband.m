function x_overlap = comm2d_reduce_wideband(x_overlap, overlap_q, Qy, Qx, Qc)
% comm2d_update_reduce: perform additive reduction of the duplicated area 
% (overlap regions, referred to as ghost cells) in a 2D image tessellation
% (2D communication grid).
%-------------------------------------------------------------------------%
%%
% Input:
%
% > x_overlap             spatial facet considered (with ghost cells)
% > overlap_q             size of the left and top facet extensions 
%                         (ghost cells) [2 ,1] 
% > Qy                    number of spatial facets along the y axis [1]
% > Qx                    number of spatial facets along the x axis [1]
% > Qc                    number of spectral facets, along the z axis [1]
%
% Output:
%
% < x_overlap             updated spatial facet
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% [22/03/2019] final debug, code ok
% [03/04/2019] add basic support for wideband, see if further modifications 
% are needed later on (for the spectral splitting)
%-------------------------------------------------------------------------%
%%

% communications to aggregate information
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
src_diag = [];

% data to send
data_vert = [];
data_horz = [];
data_diag = [];

% define communications (to N, W, NW)
if qy > 1
    % N communication (i, (qy-1, qx)) -> (i, q = (qx-1)*Qy + qy-1)
    dest_vert = get_index(i, (qx-1)*Qy + qy-1);
    data_vert = x_overlap(1:overlap_q(1), 1:end, :);
    if qx > 1
        % NW communication (qy-1, qx-1) -> q = (qx-2)*Qy + qy-1
        dest_diag = get_index(i, (qx-2)*Qy + qy-1);
        data_diag = x_overlap(1:overlap_q(1), 1:overlap_q(2), :); 
        % another set of variables overlap_q will be needed for the update of the ghost cells (adjoint communicationss)
    end
end

if qx > 1
    % W communication (qy, qx-1) -> q = (qx-2)*Qy + qy
    dest_horz = get_index(i, (qx-2)*Qy + qy);
    data_horz = x_overlap(1:end, 1:overlap_q(2), :);
end

% define receptions (from S, E, SE)
if qy < Qy
    % S reception (qy+1, qx) -> q = (qx-1)*Qy + qy+1
    src_vert = get_index(i, (qx-1)*Qy + qy+1);
    if qx < Qx
        % SE reception (qy+1, qx+1) -> q = qx*Qy + qy+1
        src_diag = get_index(i, qx*Qy + qy+1);
    end
end

if qx < Qx
    % E reception (qy, qx+1) -> q = qx*Qy + qy
    src_horz = get_index(i, qx*Qy + qy);
end

% is there a way to do everything at once? (i.e., reduce the
% synchronization phasea induced by the three calls to labSendReceive?)
% see if this can be as efficient as the gop operation (to be confirmed)
% vertical communications
rcv_vert_data = labSendReceive(dest_vert, src_vert, data_vert);
% horizontal communications
rcv_horz_data = labSendReceive(dest_horz, src_horz, data_horz);
% diagonal communications
rcv_diag_data = labSendReceive(dest_diag, src_diag, data_diag);

% update portions of the overlapping facet with the received data (aggregate and sum)
if ~isempty(rcv_vert_data) % from S
    x_overlap(end-size(rcv_vert_data, 1)+1:end, 1:end, :) = ...
        x_overlap(end-size(rcv_vert_data, 1)+1:end, 1:end, :) + rcv_vert_data;
end

if ~isempty(rcv_horz_data) % from E
    x_overlap(1:end, end-size(rcv_horz_data, 2)+1:end, :) = ...
        x_overlap(1:end, end-size(rcv_horz_data, 2)+1:end, :) + rcv_horz_data;
end

if ~isempty(rcv_diag_data) % from SE
    x_overlap(end-size(rcv_diag_data, 1)+1:end, end-size(rcv_diag_data, 2)+1:end, :) = ...
        x_overlap(end-size(rcv_diag_data, 1)+1:end, end-size(rcv_diag_data, 2)+1:end, :) + rcv_diag_data;
end

end
