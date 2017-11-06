function dy = Six(t,y)
    y = reshape(y,6,1);
    x1 = y(1);
    x2 = y(2);
    u1 = y(3);
    u2 = y(4);
    u3 = y(5);
    u4 = y(6);
    dy = zeros(6,1);
    dy(1) = x2;
    dy(2) = -(x1 + x2) * u2;
    dy(3) = (x1 + x2)*(x2 - x1 * u2 - x2 * u2) * u2;
    dy(4) = -(x1 + x2) * (x2 - x1 * u2 - x2 * u2) * u2^3;
    dy(5) = x1 * x2 * u4;
    dy(6) = -x1 * x2 * u4^3;