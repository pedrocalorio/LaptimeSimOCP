function [z,pack] = packDecVar(x,u)

nGrid = size(x,2);
nState = size(x,1);
nControl = size(u,1);

xCol = reshape(x, nState*nGrid, 1);
uCol = reshape(u, nControl*nGrid, 1);

indz = reshape(1:numel(u)+numel(x),nState+nControl,nGrid);

% index of state and control variables in the decVar vector
xIdx = indz(1:nState,:);
uIdx = indz(nState+(1:nControl),:);

% decision variables are indexed so that the defects gradients appear as a banded matrix
z = zeros(numel(indz),1);
z(xIdx(:),1) = xCol;
z(uIdx(:),1) = uCol;

pack.nGrid = nGrid;
pack.nState = nState;
pack.nControl = nControl;
pack.xIdx = xIdx;
pack.uIdx = uIdx;

end