constants_rimless; 

%% Symbolics
syms q1 q2 q1_dot q2_dot

syms g ph y p T

q = [q1 ; q2 ] ;

q_dot = [q1_dot; q2_dot];

x = [q; q_dot];

%% Create DH
q = sym('q', [n,1]);
dh = sym(DH);

dh(:,end) = DH(:,end) + q;

%% Transformation
dh_ = sym('dh', [1,4]);

a_i = dh_(2);
alpha_i = dh_(1);
d_i = dh_(3);
q_i = dh_(4);

T_i= [cos(q_i) -sin(q_i)*cos(alpha_i) sin(q_i)*sin(alpha_i) a_i*cos(q_i);...
      sin(q_i) cos(alpha_i)*cos(q_i) -cos(q_i)*sin(alpha_i) a_i*sin(q_i);...
      0 sin(alpha_i) cos(alpha_i) d_i;...
      0 0 0 1];

matlabFunction(T_i, 'File', 'symbolic_functions/Create_Transformation', 'Vars', {dh_} );

%% Transformation


for i=1:n
    T_local(:,:,i) = Create_Transformation(dh(i,:));
end

T_world(:,:,1)  = T_local(:,:,1);


for i=2:(n)
    
    T_world(:,:,i)= T_world(:,:,i-1)*T_local(:,:,i);
    
end

matlabFunction(T_world, 'File', 'symbolic_functions/get_0Tn', 'Vars', {q});

%% Postion 

for i=1:n
P(1:3,i)=T_world(1:3,4,i);
end

%% Jacobian



for i=1:n
    
    x1=P(1,i);
    y1=P(2,i);
    J_n_linear(1:2,1:n,i)= jacobian([x1;y1],[q.']) ;
    
end
J_n_linear(n+1,:)=0;

matlabFunction(J_n_linear, 'File', 'symbolic_functions/get_J_lin', 'Vars', {q});


for i=1:n
    if i>1
        J_n_w(1:3,i-1,i) = J_n_w(1:3,i-1,i-1);
    end
    J_n_w(1:3,i,i)   = (T_world(1:3,3,i));
end

matlabFunction(J_n_w, 'File', 'symbolic_functions/get_J_ang', 'Vars', {q});


for i=1:n
    J_n(:,:,i) = [J_n_linear(:,:,i) ; J_n_w(:,:,i)];
end

matlabFunction(J_n, 'File', 'symbolic_functions/get_J', 'Vars', {q});



%% Linear KE

syms m1 m2 D_linear  K_E_linear a D_q 



    D_linear = 0;

for i=1:n
    
    a = m(i)*(J_n_linear(:,:,i)')*J_n_linear(:,:,i);
    
    D_linear = D_linear + a;
    
end
    % Trig substitutions
    D_linear = subs( D_linear,[conj(q1),conj(q2),conj(l1),conj(l2)],[q1,q2,l1,l2]);
    D_linear = subs( D_linear,[sin(q1)*sin(q1) + cos(q1)*cos(q1),sin(q2)*sin(q2) + cos(q2)*cos(q2),sin(q1)*cos(q2) + cos(q1)*sin(q2),sin(q1)*cos(q2) - cos(q1)*sin(q2),sin(q2)*cos(q1) - cos(q2)*sin(q1),cos(q1)*cos(q2) - sin(q1)*sin(q2),cos(q1)*cos(q2) + sin(q1)*sin(q2),cos(q2)*cos(q1) - sin(q2)*sin(q1)],[1,1,sin(q1 + q2),sin(q1 - q2),sin(q2 + q1),cos(q1 + q2),cos(q1 - q2),cos(q1 + q2)]);        %%trigonometric substitutions


% K_E_linear = (1/2) * q_dot' * D_linear * q_dot;

%% Angular 
  

syms D_angular K_E_angular b

    K_w = sym('K_w', [1 2]);
    D_angular = 0;
    D_q = D_linear + D_angular;

for i=1:2
    
    
    b =  K_w(i);
    b = (J_n_w(:,:,i)' *  (T_world(1:3,1:3,i)) * I1 *  T_world(1:3,1:3,i)' * J_n_w(:,:,i));
    
    D_angular = D_angular + b;
end
    % Trig substitutions
    D_angular = subs( D_angular,[conj(q1),conj(q2),conj(l1),conj(l2)],[q1,q2,l1,l2]);
    D_angular = subs( D_angular,[sin(q1)*sin(q1) + cos(q1)*cos(q1),sin(q2)*sin(q2) + cos(q2)*cos(q2),sin(q1)*cos(q2) + cos(q1)*sin(q2),sin(q1)*cos(q2) - cos(q1)*sin(q2),sin(q2)*cos(q1) - cos(q2)*sin(q1),cos(q1)*cos(q2) - sin(q1)*sin(q2),cos(q1)*cos(q2) + sin(q1)*sin(q2),cos(q2)*cos(q1) - sin(q2)*sin(q1)],[1,1,sin(q1 + q2),sin(q1 - q2),sin(q2 + q1),cos(q1 + q2),cos(q1 - q2),cos(q1 + q2)]);        %%trigonometric substitutions

%     K_E_angular = (1/2) * q_dot' * D_angular * q_dot;

%% D matrix

    D_q = D_linear + D_angular;
    
    matlabFunction(D_q, 'File', 'symbolic_functions/get_D', 'Vars', {q});

%% Total KE

    K_E = (1/2) * q_dot' * D_q * q_dot  ;
    
    matlabFunction(K_E, 'File', 'symbolic_functions/get_KE', 'Vars', {x});
%% Corriolis
c(1:3,1:3)=sym(zeros(3));
for i=1:n
    for j=1:n
        for k=1:n
            c(i,j,k)= (1/2) * (c(i,j,k)+ (jacobian( D_q(k,j) , q(i)) + (jacobian( D_q(k,i) , q(j)) - (jacobian( D_q(i,j) , q(k) )))));
        end
    end
    
end

matlabFunction(c, 'File', 'symbolic_functions/get_C', 'Vars', {q});

%% Potential Energy
Y_com=sym('ycom',[1 2]);
Q=0;
for i=1:n
    y = Y_com(i);
    Q = Q+q(i);
    Y_com(i)= L_com(i)*sin(Q);
    
end

P_i = sym('P', [1 2]);
phi_i = sym('phi', [1 2]);
Pi=0;
for i=1:n
    p=P_i(i);
    
    p= m(i)*grav*Y_com(i);
    Pi = Pi + p;
    
end

phi_q = 0;

for i=1:n
    ph = phi_i(i);
    ph = jacobian(Pi,q(i));
    phi_q = phi_q + ph;
end
matlabFunction(phi_q, 'File', 'symbolic_functions/get_Phi', 'Vars', {q});


