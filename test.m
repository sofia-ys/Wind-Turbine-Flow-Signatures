% basic arithmetic
a = 5;  % if you don't use a semicolon, the output will be displayed
b = 3;  % without a semicolon with would say "b = 3" in the command window
sum = a + b  % since this does not have a semicolon the full statement is shown in the command window
difference = a - b
product = a * b
quotient = a / b 

% vectors and matrices
v = [1, 2, 3, 4, 5]  % row vector
w = [1; 2; 3; 4; 5; 6]  % column vector
A = [1 2 3; 4 5 6; 7 8 9]  % 3x3 matrix
B = A'  % this computes the transpose of A

% plotting a function
x = linspace(0, 10, 100);  % creates 100 points between 0 and 10
y = sin(x);
plot(x,y)
xlabel('x')
ylabel('y')
title('Sine Wave')
grid on

% loops and conditions
for i = 1:5  % for i in range 1 to 5 with a step of 1: [1, 2, 3, 4, 5]
    fprintf('Iteration %d\n', i);  % fprintf is formatted printing, %d specifies a decimal integer
end  % always have to include an end statement

x = 10;
if x > 5
    disp('x is greater than 5')  % print without formatting
else
    disp('x is 5 or less')
end

% defining and calling a function
function y = squareNumber(x)  % function called squareNumber takes one input x and returns one output y
    y = x^2;  % without a semicolon this would also print y = 16
end

result = squareNumber(4)

function result = addMultiply(a, b, c)  % function with multiple inputs
    result = (a + b) * c;
end

x = addMultiply(2, 3, 4)  % (2 + 3) * 4 = 20

function [sumResult, productResult] = sumAndProduct(a,b)  % function with multiple inputs and outputs
    sumResult = a + b;
    productResult = a * b;
end

[sumVal, prodVal] = sumAndProduct(4, 5)  % [] are used to capture multiple outputs