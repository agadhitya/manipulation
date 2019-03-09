function J_n = get_J(in1)
%GET_J
%    J_N = GET_J(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    08-Mar-2019 19:15:20

q1 = in1(1,:);
q2 = in1(2,:);
t2 = sin(q1);
t3 = cos(q1);
t4 = sin(q2);
t5 = cos(q2);
t6 = t3.*t5;
J_n = reshape([-t2,t3,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,-t2-t2.*t5-t3.*t4,t3+t6-t2.*t4,0.0,0.0,0.0,1.0,-t2.*t5-t3.*t4,t6-t2.*t4,0.0,0.0,0.0,1.0],[6,2,2]);