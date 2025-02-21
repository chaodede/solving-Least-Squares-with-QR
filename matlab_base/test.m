% V1 = V0 - pinv(J) * e; update unknowns
% change equantion to Least-Squares Problem
% b =J * V0 - e 
% J * V1 = b
J = [-86.061506900000000 -0.020595972000000 0 0; 
    -255.507089802444 -3.53229983209699 -635.836267038566 5.72169364259395E-21;
    -253.17974913249 -2.51170386044137 -2290.73769949779 -10.0998528315302;
    4.28903178185595E-24 4.25499184707931E-26 -1653.90143245922 -11.0998528315302
    ];
b = [2.0492325907922506;
    -54188.500682648504;
    -188565.40871117834;
    -134700.43885481835];
V1 = [-0.1505;
    529.3641;
    82.3237;
    -134.0421];
R = J;
Q = eye(4);
P = eye(4);
% column 1
lens = [R(:,1)'*R(:,1) R(:,2)'*R(:,2) R(:,3)'*R(:,3) R(:,4)'*R(:,4)];
disp(lens);
R = R(:,[3 2 1 4]);
P = P(:,[3 2 1 4]);
lens = lens*P;
v = get_house(R(:,1),1,4);
R = R-v*(v'*R);
Q = Q - (Q*v)*v';
disp(R);
% column 2
lens(2:4) = lens(2:4) - [R(1,2)^2 R(1,3)^2 R(1,4)^2];
disp(lens);
R = R(:,[1 3 2 4]);
P = P(:,[1 3 2 4]);
lens = lens*P;
v = get_house(R(:,2),2,4);
R = R-v*(v'*R);
Q = Q - (Q*v)*v';
disp(R);
% column 3
lens(3:4) = lens(3:4) - [R(2,3)^2 R(2,4)^2];
disp(lens);
v = get_house(R(:,3),3,4);
R = R-v*(v'*R);
Q = Q - (Q*v)*v';
disp(R);
% if there is square matrix, no need to do the last column
% column 4 
% v = get_house(R(:,4),4,4);
% R = R-v*(v'*R);
% Q = Q - (Q*v)*v';
% disp(R);

% We are done and we can check with
check_diff = J*P - Q*R;
Ri = inv(R);
check_diff3 = Ri*Q';
z = upper_solve(R,(Q'*b));
x = P*z;
disp(x);

% If we have the right solution then A*x and b - A*x should be orthogonal
check_diff2 = (J*x)'*(b-J*x);



