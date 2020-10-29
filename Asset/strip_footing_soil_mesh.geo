SetFactory("OpenCASCADE");

// external rectangle
Point(1) = { -5.0,  0.25, 0.0, 0.125 };
Point(2) = {  5.0,  0.25, 0.0, 0.125 };
Point(3) = { -5.0,  0.0, 0.0, 0.125 };
Point(4) = {  5.0,  0.0, 0.0, 0.125 };
Point(5) = { -5.0, -3.0, 0.0, 0.125 };
Point(6) = {  5.0, -3.0, 0.0, 0.125 };
Point(7)  = { -5.0, -3.5, 0.0, 0.3 };
Point(8)  = {  5.0, -3.5, 0.0, 0.3 };
Point(9)  = { -5.0, -5.0, 0.0, 0.3 };
Point(10) = {  5.0, -5.0, 0.0, 0.3 };

Line(1) = { 1, 2 };
Line(2) = { 1, 3 };
Line(3) = { 2, 4 };
Line(4) = { 3, 4 };
Line(5) = { 3, 5 };
Line(6) = { 4, 6 };
Line(7) = { 5, 6 };
Line(8) = { 5, 7 };
Line(9) = { 6, 8 };
Line(10) = { 7, 8 };
Line(11) = { 7, 9 };
Line(12) = { 8, 10 };
Line(13) = { 9, 10 };

Curve Loop(1) = { 1,  3,  -4, -2  };
Curve Loop(2) = { 4,  6,  -7, -5  };
Curve Loop(3) = { 7,  9,  -10, -8  };
Curve Loop(4) = { 10, 12, -13, -11  };

Plane Surface(1) = { 1 };
Plane Surface(2) = { 2 };
Plane Surface(3) = { 3 };
Plane Surface(4) = { 4 };
