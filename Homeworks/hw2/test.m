
syms s

A = [0 1; 0 -1];
B = [ 0; 1];
C = [10 5];

TF1 = C*inv(s*eye(size(A))-A)*B

A = [0 1; 0 -1];
B = [0; 5];
C = [2 1];

TF2 = C*inv(s*eye(size(A))-A)*B