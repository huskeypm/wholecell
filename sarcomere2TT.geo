ll = 0.2;
lc = 0.1;
tRad = 0.25;
Point(1) = {0.0, 1.5, 0, lc};  
Point(2) = {2.0, 1.5, 0, lc};  
Point(3) = {2.0, 0.0, 0, lc};  
Point(4) = {0.0, 0.0, 0, lc};  
Point(5) = {0.0      , 0.75 - tRad, 0, lc};  
Point(6) = {tRad     , 0.75       , 0, lc};  
Point(16) ={0.00     , 0.75       , 0, lc};  
Point(7) = {0.0      , 0.75 + tRad, 0, lc};  
Point(8) = {2.0      , 0.75 - tRad, 0, lc};  
Point(9) = {2. - tRad, 0.75       , 0, lc};  
Point(19) ={2.0      , 0.75       , 0, lc};  
Point(10) ={2.0      , 0.75 + tRad, 0, lc};  

Line(1) = {1, 2};
Line(2) = {10, 2};
Line(3) = {8, 3};
Line(4) = {3, 4};
Line(5) = {4, 5};
Line(6) = {7, 1};
Circle(7) = {7, 16, 6};
Circle(8) = {6, 16, 5};
Circle(9) = {10, 19, 9};
Circle(10) = {9, 19, 8};

/* define surface */
Line Loop(11) = {6, 1, -2, 9, 10, 3, 4, 5, -8, -7};

/* extrude surface */
Plane Surface(12) = {11};
Extrude {0, 0, 1} {
  Surface{12};
}
