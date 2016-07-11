clc;

%-- initialization and data input part --%
% set total number of sample
N = size(ax,1);
% set sampling period
dt = 0.01;
% set sampling time
T = 0:dt:N*dt-dt;
% set gravity
g = 9.81;
% input measured acceleration, angular rate, attitude
measuredA = [ax,ay,az]' * g;
measuredW = [wx,wy,wz]';
measuredAngle = [roll,pitch,yaw]';

%--  trajectory reconstruction with normal kalman filter and MAUKF --%

% initialize state vector, distance estimation, yaw angle and error of yaw angle
% OS = (wx,wy,wz,q0,q1,q2,q3)
OS = zeros(7,N);
% PS = (p1x, v1x, a1x, p1y, v1y, a1y, p1z, v1z, a1z?
PS = zeros(9,N);
% set time constant for body angular velocity
tao1 = 0.5;
tao2 = 0.5;
tao3 = 0.5;
D1 = 0.4;
D2 = 0.4;
D3 = 0.4;
% set system input of orientation estimation
OB = 0;
Ou = [0,0,0,0,0,0,0]';
% set orientation measurement matrix
OC = eye(7);
% set orientation process noise covariance
OQ = diag([D1/2/tao1*(1-exp(-2*dt/tao1)), D2/2/tao2*(1-exp(-2*dt/tao2)), D3/2/tao3*(1-exp(-2*dt/tao3)), 0, 0, 0, 0]);
% set orientationmeasure noise covariance
OR = diag([0.01,0.01,0.01,0.0001,0.0001,0.0001,0.0001]);
% set orientation forcast error covariance
OP = 0.1 * eye(7);
% set system input of position estimation
Pu = zeros(9,1);
PB = 0;
% set position process noise(m/s^2)
w = 0.02;
% set position measure noise(m/s^2)
z = 0.1;
% set position process noise covariance
PQ = w^2 .* diag([dt*dt,dt,1,dt*dt,dt,1,dt*dt,dt,1]);
% set position measure noise covariance
PR = z^2 .* eye(3);
% set position forcast error covariance
PP = 0.5 * diag([dt*dt,dt,1,dt*dt,dt,1,dt*dt,dt,1]);
% initialize stop count
stopCount = zeros(3,1);
for i = 1 : N
    if i == 1      
        OS(:,i) = [measuredW(1:3,i); EulerAngleToQuaternion(measuredAngle(1:3,i))];
    else
        % calculate orientation state transition matrix
        OA = calculateOrientationTransitionMatrix(OS(:,i-1), tao1, tao2, tao3, dt);
        % set orientation state measurement matrix
        OC = eye(7);
        % do Kalman filter for orientation
        [OS(:,i),OP] = kalmanFilter(OS(:,i-1), OP, OA, OB, Ou, OC, [measuredW(1:3,i);EulerAngleToQuaternion(measuredAngle(1:3,i))], OQ, OR);
        % renormalize quaternion
        OS(4:7,i) = OS(4:7,i)/(sqrt(sum(OS(4:7,i) .* OS(4:7,i))));
    end
    % calculate rotation matrix from quaternion
    RM = QuaternionToRotationMatrix(OS(4:7,i));
    if i == 1
        PS(:,i) = [0, 0, measuredA(1,i), 0, 0, measuredA(2,i), 0, 0, measuredA(3,i)];
    else
        % calculate position state transiton matrix
        PA = [1,dt,RM(1,1)*dt*dt/2,0,0,RM(1,2)*dt*dt/2,0,0,RM(1,3)*dt*dt/2; 0,1,RM(1,1)*dt,0,0,RM(1,2)*dt,0,0,RM(1,3)*dt; 0,0,1,0,0,0,0,0,0; ...
              0,0,RM(2,1)*dt*dt/2,1,dt,RM(2,2)*dt*dt/2,0,0,RM(2,3)*dt*dt/2; 0,0,RM(2,1)*dt,0,1,RM(2,2)*dt,0,0,RM(2,3)*dt; 0,0,0,0,0,1,0,0,0; ...
              0,0,RM(3,1)*dt*dt/2,0,0,RM(3,2)*dt*dt/2,1,dt,RM(3,3)*dt*dt/2; 0,0,RM(3,1)*dt,0,0,RM(3,2)*dt,0,1,RM(3,3)*dt; 0,0,0,0,0,0,0,0,1];
        % set position measurement matrix
        PC = [0,0,1,0,0,0,0,0,0; 0,0,0,0,0,1,0,0,0; 0,0,0,0,0,0,0,0,1];
        % do normal Kalman filter for position
        [PS(:,i),PP] = kalmanFilter(PS(:,i-1), PP, PA, PB, Pu, PC, measuredA(:,i), PQ, PR);
    end
    % calculate device acceleration to earth
    
    aToEarth = RM * measuredA(1:3,i);
    
%     if i == 2500 || i == 4500 || i == 7000
%         PS(2,i) = 0;
%         PS(5,i) = 0;
%         PS(8,i) = 0;
%     end
    for j=1:size(aToEarth,1)
        if abs(aToEarth(j)) < 0.02 * g * 2
            stopCount(j) = stopCount(j)+1;
        else
            stopCount(j) = 0;
        end
        if stopCount(j) == 100
            PS(3*j-1,i) = 0;
            stopCount(j) = 0;
        end
    end
end

PS(7,1:6)
%-- consequence exhibitoin part --%
figure
plot3(PS(1,:), PS(4,:), PS(7,:),'b');
xlabel('x/m');
ylabel('y/m');