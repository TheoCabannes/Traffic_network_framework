A = [2 1 1; 1 1 0; 1 0 1]
B = [2.25 3 3]

G = [1 0.1; 2 0; 0.25 0; 2 0; 1 0.1]


nb_iter = 13
flow = zeros(nb_iter, 3)
tt = zeros(nb_iter, 3)
obj = zeros(nb_iter,1)
x = zeros(nb_iter, 1)

for tau=[1:nb_iter]
% test of social optimum and Nash equilibrium using intergral formulation
c_middle = (tau - 1) / 2
c = [11 11 c_middle 11 11]
x(tau) = c_middle
cvx_begin
    variables f(5)
	f(3) == f(1) - f(4) % the path flow constraint
    f(3) == f(5) - f(2) % the path flow constraint
    f(1)+f(2)==10 % the demand constraint
    f >= 0 % flow should be positive
    f <= c'
    % nash equilibrium
    minimize ((G(1:5, 1)' * f) + (1/2 * G(1:5,2)' * (f .* f)))
cvx_end

obj(tau) = ((G(1:5, 1)' * f) + (1/2 * G(1:5,2)' * (f .* f)))
h = [f(3) f(4) f(2)]
flow(tau, :) = h
tt(tau, :) = (A * h' / sum(h)) + B'

end

plot(x, flow)
title("The path flow as a function of Capacity")
xlabel("Lambda")
ylabel("Path flow")
figure
plot(x, tt)
title("The path travel time as a function of Capacity")
xlabel("Capacity")
ylabel("Path travel time")
figure
plot(x, obj)
title("The objective function as a function of Capacity")
xlabel("Capacity")
ylabel("Objective function")
