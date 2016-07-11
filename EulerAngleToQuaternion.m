function q = EulerAngleToQuaternion(angle)
q = [cos(angle(3)/2) * cos(angle(2)/2) * cos(angle(1)/2) + sin(angle(3)/2) * sin(angle(2)/2) * sin(angle(1)/2);
     cos(angle(3)/2) * cos(angle(2)/2) * sin(angle(1)/2) - sin(angle(3)/2) * sin(angle(2)/2) * cos(angle(1)/2);
     cos(angle(3)/2) * sin(angle(2)/2) * cos(angle(1)/2) + sin(angle(3)/2) * cos(angle(2)/2) * sin(angle(1)/2);
     sin(angle(3)/2) * cos(angle(2)/2) * cos(angle(1)/2) - cos(angle(3)/2) * sin(angle(2)/2) * sin(angle(1)/2)];
q = q/(sqrt(sum(q .* q)));
end