t = 0.125;
theta = [-75.6 -46.8 45 -34.2]; % In Degrees
theta = theta.*(pi/180); % Conversion to Radians - I my system Default is radians

E1 = 131;
E2 = 9.8;
G12 = 5.8;
nu_12 = 0.22;

S = [1/E1 -nu_12/E1 0; -nu_12/E1 1/E2 0; 0 0 1/G12];

Q_local = S^-1;

T = @(x) [(cos(x))^2 (sin(x))^2 2*cos(x)*sin(x); (sin(x))^2 (cos(x))^2 -2*cos(x)*sin(x); -cos(x)*sin(x) cos(x)*sin(x) (cos(x))^2 - (sin(x))^2];

%Q_global = (inv(T))*(Q_local)*(inv(transpose(T)));

A = zeros(3);
B = zeros(3);
D = zeros(3);

for k = 1:4
   z_k = (3-k)*0.125;
   z_k1 = z_k - 0.125;
   Q_global = ((inv(T(theta(k))))*(Q_local)*(inv(transpose(T(theta(k))))));
   A = A + 0.125*Q_global;
   B = B + (0.125*(z_k + z_k1))*Q_global;
   D = D + (0.125*(z_k^2 + z_k1^2 + z_k*z_k1))*Q_global;
end

sol = [A B;B D]