SetFactory("OpenCASCADE");

// external rectangle
Point(1) = { 0.0, 0.0, 0.0, 0.025 };
Point(2) = { 0.2, 0.0, 0.0, 0.025 };
Point(3) = { 0.2, 1.0, 0.0, 0.025 };
Point(4) = { 0.0, 1.0, 0.0, 0.025 };
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};