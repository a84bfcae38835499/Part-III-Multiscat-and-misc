clear;

A = [1 2 3;
     4 5 6;
     7 8 9];

B = transpose(A);
B = rot90(A);
B(1,1) = 0;


disp(AreCyclicBoundaryMatriciesEqual(A,B))