SetFactory("OpenCASCADE");

// external rectangle
Point(1) = { 0.0, 0.5, 0.0, 0.05 };
Point(2) = { 3.0, 0.5, 0.0, 0.05 };
Point(3) = { 3.5, 0.5, 0.0, 0.15 };
Point(4) = { 5.0, 0.5, 0.0, 0.15 };

Point(5) = { 0.0, 0.0, 0.0, 0.05 };
Point(6) = { 3.0, 0.0, 0.0, 0.05 };
Point(7) = { 3.5, 0.0, 0.0, 0.15 };
Point(8) = { 5.0, 0.0, 0.0, 0.15 };

Point(9) = { 0.0, -3.0, 0.0, 0.05 };
Point(10) = { 3.0, -3.0, 0.0, 0.05 };

Point(11) = { 0.0, -3.5, 0.0, 0.15 };
Point(12) = { 3.5, -3.5, 0.0, 0.15 };
Point(13) = { 5.0, -3.5, 0.0, 0.15 };

Point(14) = { 0.0, -5.0, 0.0, 0.15 };
Point(15) = { 5.0, -5.0, 0.0, 0.15 };

Line(1) = { 1, 2 };
Line(2) = { 2, 3 };
Line(3) = { 3, 4 };

Line(4) = { 5, 6 };
Line(5) = { 6, 7 };
Line(6) = { 7, 8 };

Line(7) = { 9, 10 };
Line(8) = { 11, 12 };
Line(9) = { 12, 13 };
Line(10) = { 14, 15 };

Line(11) = { 1, 5 };
Line(12) = { 5, 9 };
Line(13) = { 9, 11 };
Line(14) = { 11, 14 };

Line(15) = { 2, 6 };
Line(16) = { 6, 10 };

Line(17) = { 3, 7 };
Line(18) = { 7, 12 };

Line(19) = { 4, 8 };
Line(20) = { 8, 13 };
Line(21) = { 13, 15 };

Curve Loop(1) = { 1, 15, 4, 11 };
Curve Loop(2) = { 2, 17, 5, 15 };
Curve Loop(3) = { 3, 19, 6, 17 };
Curve Loop(4) = { 4, 16, 7, 12 };
Curve Loop(5) = { 5, 18, 8, 13, 7, 16 };
Curve Loop(6) = { 6, 20, 9, 18 };

Curve Loop(7) = { 8, 9, 21, 10, 14 };

Plane Surface(1) = { 1 };
Plane Surface(2) = { 2 };
Plane Surface(3) = { 3 };
Plane Surface(4) = { 4 };
Plane Surface(5) = { 5 };
Plane Surface(6) = { 6 };
Plane Surface(7) = { 7 };