
% x = zeros(nState,nGrid);
% for i=1:nState
% %    x(i,:) = zeros(1,nState); 
%    for j=1:nGrid
%        x(i,j) = z(1,i:i+nGrid);
%    end
% 
% end


vx = z(1,(3-1)*nGrid+1:3*nGrid);
vy = z(1,(4-1)*nGrid+1:4*nGrid);

x = [vx;vy];


