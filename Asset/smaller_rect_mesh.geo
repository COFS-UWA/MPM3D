SetFactory("OpenCASCADE");

// external rectangle
Point(1) = {  0.0, 0.0,  0.0, 1.0 };
Point(2) = { 20.0, 0.0,  0.0, 1.0 };
Point(3) = { 20.0, 15.0, 0.0, 1.0 };
Point(4) = {  0.0, 15.0, 0.0, 1.0 };
Point(5) = {  0.0, 20.0, 0.0, 1.0 };
Point(6) = { 20.0, 20.0, 0.0, 1.0 };
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {4, 5};
Line(6) = {5, 6};
Line(7) = {6, 3};
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {3, 5, 6, 7};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

//field[50] = matheval; //generate field
//field[50].f = "5.0";
//background field = 50;