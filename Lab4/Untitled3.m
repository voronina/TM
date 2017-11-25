clear;

fid = fopen('results.txt');
input = cell2mat(textscan(fid,'%f'));
fclose(fid);
k = 1;
for i = 1:2:length(input)
    input_x(k) = input(i);
    input_y(k) = input(i + 1);
    k = k + 1;
end

% Границы графика
LEFT = -10;
RIGHT = 10;

% Исходные параметры
A1 = 0.2;
A2 = 1;
A3 = -0.1;

% Задание исходной функций
fun = @(a) and( a >= -8, a < -4 ).*( A1 ) ...
    + and( a >= -4, a < 1 ).*( A2 ) ...
    + and( a >= 1, a <= 9 ).*( A3 );

exp = 10;
length = 2^exp;
step = abs(RIGHT - LEFT)/length;

range = LEFT:step:RIGHT-step;
points = fun(range);
Y = fft(fun(range));

plot(input_x, input_y, input_x,abs(fftshift(Y)/length),'--');
legend('FILE','MATLAB');
grid on;

clear A1 A2 A3 fid fun i input ans k;
