function OA = calculateOrientationTransitionMatrix(OS,tao1,tao2,tao3,dt)
% OA = [diag([exp(-dt/tao1), exp(-dt/tao2), exp(-dt/tao3)]), zeros(3,4); ...
%       -OS(5)*dt/2, -OS(6)*dt/2, -OS(7)*dt/2, 1, -OS(1)*dt/2, -OS(2)*dt/2,-OS(3)*dt/2; ...
%       OS(4)*dt/2, -OS(7)*dt/2, OS(6)*dt/2, OS(1)*dt/2, 1, OS(3)*dt/2, -OS(2)*dt/2; ...
%       OS(7)*dt/2, OS(4)*dt/2, -OS(5)*dt/2, OS(2)*dt/2, -OS(3)*dt/2, 1, OS(1)*dt/2; ...
%       -OS(6)*dt/2, OS(5)*dt/2, OS(4)*dt/2, OS(3)*dt/2, OS(2)*dt/2, -OS(1)*dt/2, 1];
OA = [diag([exp(-dt/tao1), exp(-dt/tao2), exp(-dt/tao3)]), zeros(3,4); ...
      zeros(4,3),[1, -OS(1)*dt/2, -OS(2)*dt/2,-OS(3)*dt/2; OS(1)*dt/2, 1, OS(3)*dt/2, -OS(2)*dt/2; ...
      OS(2)*dt/2, -OS(3)*dt/2, 1, OS(1)*dt/2; OS(3)*dt/2, OS(2)*dt/2, -OS(1)*dt/2, 1]];
end