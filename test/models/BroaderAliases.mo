model BroaderAliases
Real x1(start = 0);
Real x2(start = 0);
Real z(start = 15, max = 25, min = 10);
Real f(start = 5, max = 10, min = 5);
input Real u(min = -5, max = 10);
input Real u1(min = -5, max = 100);
Real k;
Real l;
equation
der(x2) = x1 + u1 + l;
der(x1) = u + z;
f = 0.2;
sin(z) - 2*f = 0;
2/k = 5*f;
l - 6*x1 = 0;
end BroaderAliases;