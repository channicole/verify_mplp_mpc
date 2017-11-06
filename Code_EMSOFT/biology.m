function dx = biology(t,x)
dx = zeros(7,1);

dx(1) =  -0.4*x(1) + 50*x(3)*x(4);
dx(2) =  0.4*x(1) - x(2);
dx(3) =  x(2) - 50*x(3)*x(4);
dx(4) =  50*x(5)*x(6) - 50*x(3)*x(4);
dx(5) =  -50*x(5)*x(6) + 50*x(3)*x(4);
dx(6) = 0.5*x(7) - 50*x(5)*x(6);
dx(7) = -0.5*x(7) + 50*x(5)*x(6);