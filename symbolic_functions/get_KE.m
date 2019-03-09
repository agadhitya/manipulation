function K_E = get_KE(in1)
%GET_KE
%    K_E = GET_KE(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    08-Mar-2019 19:15:21

q1 = in1(1,:);
q2 = in1(2,:);
q1_dot = in1(3,:);
q2_dot = in1(4,:);
t2 = q1+q2;
t3 = sin(t2);
t4 = cos(t2);
t5 = conj(q1_dot);
t6 = cos(q1);
t7 = t4+t6;
t8 = sin(q1);
t9 = t3+t8;
t10 = conj(q2_dot);
t11 = t3.*t9;
t12 = t4.*t7;
t13 = t11+t12+1.0;
K_E = q2_dot.*((t5.*t13)./2.0+(t10.*(t3.^2+t4.^2+1.0))./2.0)+q1_dot.*((t10.*t13)./2.0+(t5.*(t7.^2+t9.^2+3.0))./2.0);