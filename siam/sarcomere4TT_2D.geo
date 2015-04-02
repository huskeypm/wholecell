lc = 0.1;
tRad = 0.25;

Point(1) = {0.0, 1.5, 0, lc};  
Point(11) = {0.0+tRad, 1.5, 0, lc};  
Point(12) = {0.0, 1.5-tRad, 0, lc};  
Point(2) = {2.0, 1.5, 0, lc};  
Point(21) = {2.0-tRad, 1.5, 0, lc};  
Point(22) = {2.0, 1.5-tRad, 0, lc};  
Point(3) = {2.0, 0.0, 0, lc};  
Point(31) = {2.0-tRad, 0.0, 0, lc};  
Point(32) = {2.0, 0.0+tRad, 0, lc};  
Point(4) = {0.0, 0.0, 0, lc};  
Point(41) = {0.0+tRad, 0.0, 0, lc};  
Point(42) = {0.0, 0.0+tRad, 0, lc};  

Line(1) = {12, 42};
Line(2) = {41, 31};
Line(3) = {32, 22};
Line(4) = {21, 11};
Circle(5) = {12, 1, 11};
Circle(6) = {21, 2, 22};
Circle(7) = {32, 3, 31};
Circle(8) = {41, 4, 42};
Line Loop(9) = {4, -5, 1, -8, 2, -7, 3, -6};

Plane Surface(10) = {9};
