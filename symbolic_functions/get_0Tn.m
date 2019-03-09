function T_world = get_0Tn(in1)
%GET_0TN
%    T_WORLD = GET_0TN(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    08-Mar-2019 19:15:20

q1 = in1(1,:);
q2 = in1(2,:);
t2 = cos(q1);
t3 = sin(q1);
t4 = sin(q2);
t5 = cos(q2);
t6 = t2.*t5;
t7 = t6-t3.*t4;
t8 = t2.*t4;
t9 = t3.*t5;
T_world = reshape([t2,t3,0.0,0.0,-t3,t2,0.0,0.0,0.0,0.0,1.0,0.0,t2,t3,0.0,1.0,t7,t8+t9,0.0,0.0,-t2.*t4-t3.*t5,t7,0.0,0.0,0.0,0.0,1.0,0.0,t2+t6-t3.*t4,t3+t8+t9,0.0,1.0],[4,4,2]);