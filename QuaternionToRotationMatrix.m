function C = QuaternionToRotationMatrix(q)
    C(1,1) = q(1) * q(1) + q(2) * q(2) - q(3) * q(3) - q(4) * q(4);
    C(1,2) = 2 * (q(2) * q(3) - q(1) * q(4));
    C(1,3) = 2 * (q(2) * q(4) + q(1) * q(3));
    C(2,1) = 2 * (q(2) * q(3) + q(1) * q(4));
    C(2,2) = q(1) * q(1) - q(2) * q(2) + q(3) * q(3) - q(4) * q(4);
    C(2,3) = 2 * (q(3) * q(4) - q(1) * q(2));
    C(3,1) = 2 * (q(2) * q(4) - q(1) * q(3));
    C(3,2) = 2 * (q(3) * q(4) + q(1) * q(2));
    C(3,3) = q(1) * q(1) - q(2) * q(2) - q(3) * q(3) + q(4) * q(4);
end
