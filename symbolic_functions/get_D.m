function D_q = get_D(in1)
%GET_D
%    D_Q = GET_D(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    08-Mar-2019 19:15:21

q1 = in1(1,:);
q2 = in1(2,:);
t3 = q1+q2;
t7 = cos(t3);
t8 = cos(q1);
t2 = t7+t8;
t5 = sin(t3);
t6 = sin(q1);
t4 = t5+t6;
t9 = t4.*t5;
t10 = t2.*t7;
t11 = t9+t10+1.0;
D_q = reshape([t2.^2+t4.^2+3.0,t11,t11,t5.^2+t7.^2+1.0],[2,2]);
