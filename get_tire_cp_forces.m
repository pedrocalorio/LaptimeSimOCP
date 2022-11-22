function [Fx,Fy,Mz,Cx,Cy] = get_tire_cp_forces(mfeval_output,side)

Fx = mfeval_output(:,1)';
Fy = side.*mfeval_output(:,2)';
Mz = mfeval_output(:,6)';
Cx = mfeval_output(:,30)';
Cy = -mfeval_output(:,29)';

end