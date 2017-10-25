n = 5;
m = 90;
Before = s.AllSimus{m}{1}(:, 12, n);
[shift, After] = minimize_distance(s.AllSimus{m}{1}(:, 12, n), left2);
Average = s.AvSimus{m}{1}(:, 12);



for n = 1:10
    [shift, AfterAll(n, :)] = minimize_distance(s.AllSimus{m}{1}(:, 12, n), left2);
end
figure;
hold on;
% plot(Before);
% plot(After);
plot(Average);
plot(mean(AfterAll));
legend('Average', 'AfterAll');%'Before', 'After', 