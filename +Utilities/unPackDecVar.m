function [x,u] = unPackDecVar(z,pack)

nGrid = pack.nGrid;
nState = pack.nState;
nControl = pack.nControl;

x = z(pack.xIdx);
u = z(pack.uIdx);

% make sure x and u are returned as vectors, [nState,nTime] and [nControl,nTime]
x = reshape(x,nState,nGrid);
u = reshape(u,nControl,nGrid);

end