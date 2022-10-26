function psi_out = normalize_psi(psi)

%     # use modulo operator to remove multiples of 2*pi
psi_out = sign(psi) .* mod(abs(psi), 2 * pi);


for i=1:length(psi_out)

    if psi_out(i) >= pi
        psi_out(i) = psi_out(i) - 2 * pi;
    end

    if psi_out(i) < -pi
        psi_out(i) = psi_out(i) + 2 * pi;
    end

%     psi_out(psi_out >= pi) = psi_out - 2 * pi;
%     psi_out(psi_out < -pi) = psi_out + 2 * pi;

end

% for i=2:length(psi_out)
%  difference = psi_out(i)-psi_out(i-1);
%  if difference > pi
%     psi_out(i:end) = psi_out(i:end) - 2*pi;
%  elseif difference < -pi
%     psi_out(i:end) = psi_out(i:end) + 2*pi; 
%  end
% end






