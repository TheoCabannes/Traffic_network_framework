A = [2 1 1; 1 1 0; 1 0 1]
B = [2.25 3 3]

G = [1 0.1; 2 0; 0.25 0; 2 0; 1 0.1]

nb_iter = 13
flow = zeros(nb_iter, 3)
tt = zeros(nb_iter, 3)
x = zeros(nb_iter, 1)
y = zeros(nb_iter, 1)
for tau=[1:nb_iter]
lambda = (nb_iter - tau) / 2
x(tau) = lambda
cvx_begin
    variables f(5)
	f(3) == f(1) - f(4) % the path flow constraint
    f(3) == f(5) - f(2) % the path flow constraint
    f(1)+f(2)==10 % the demand constraint
    f >= 0 % flow should be positive
    f(3) <= lambda
    % user equilibrium
    minimize ((G(1:5, 1)' * f) + (1/2 * G(1:5,2)' * (f .* f)))
cvx_end
y(tau) = (G(1:5, 1)' * f) + (1/2 * G(1:5,2)' * (f .* f))
h = [f(3) f(4) f(2)]
flow(tau, :) = h
tt(tau, :) = (A * h' / sum(h)) + B'
end

plot(x, flow)
set(gca, 'XDir','reverse')
title("The path flow as a function of capacity")
xlabel("Capacity")
ylabel("Path flow")
figure
plot(x, tt)
set(gca, 'XDir','reverse')
title("The path travel time as a function of capacity")
xlabel("Capacity")
ylabel("Path travel time")
figure
plot(x, y)
set(gca, 'XDir','reverse')
title("The objective function as a function of the capacity")
xlabel("Capacity")
ylabel("Objective function")

