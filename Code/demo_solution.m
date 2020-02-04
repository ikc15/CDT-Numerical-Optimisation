
%% *Demo: Exercise 1*

clear all; close all

%% 
% a) Rosenbrock function
f = @(x,y) 100.*(y - x.^2).^2 + (1 - x).^2;

%%
% b) Discretisation of the unit square. It can be scaled with factor $\alpha$ using simple multiplication $\alpha X$, $\alpha Y$
% Plot the function of two variables
n = 300;
x = linspace(-1,1,n+1);
h = 2/n;
y = x;
[X,Y] = meshgrid(x,y);

alpha = 2;
%figure, surf(alpha*X, alpha*Y, f(alpha*X,alpha*Y), 'EdgeColor', 'none')
%figure, surfc(alpha*X, alpha*Y, f(alpha*X,alpha*Y), 'EdgeColor', 'none')
%figure, contour(alpha*X, alpha*Y, log(f(alpha*X,alpha*Y)))
figure, contourf(alpha*X, alpha*Y, log(f(alpha*X,alpha*Y)))

%%
% c) Gradient
%
% $$df/dx = -400 (y - x^2) x - 2 (1 - x)$$
%
% $$df/dy = 200 (y - x^2)$$

%%
% c) Hessian
%
% $$d^2f/dx^2 = -400 (y - 3x^2) + 2$$
%
% $$d^2f/dxdy = -400 x $$
%
% $$d^2f/dydx = -400 x$$
%
% $$df^2/dy^2 = 200 $$

%% 
% d) Find the minimiser $x^*$ of the function $f$. Show that $x^*$ is unique and that $\nabla^2 f(x^*)$ is positive definite.

%%
% We observe that $f(x, y) \geq 0$ because is the sum of two non-negative terms.
% Also, 
%
% $$f(x, y) = 0 \Leftrightarrow y = x^2\,\textnormal{and}\,x = 1$$
%
% This is only possible if $x = 1$ and $y = 1$.
% Hence, we have $x^* = (1, 1)^T$.

%%
% The Hessian at $x^*$ is 
%
% $$ \nabla^2 f(x^*) = \left(\begin{array}{cc} 802 & -400 \\ -400 & 200 \end{array}\right)$$
%
% Hence, to verify that it is positive definite, for all $(x, y)
% \neq 0$ we must have:
% 
% $$ \begin{array}{cc}(x & y)\end{array}\left(\begin{array}{cc} 802 & -400 \\ -400 & 200 \end{array}\right)\left(\begin{array}{c} x \\ y\end{array}\right) > 0$$
%
% This is equivalent to
%
% $$ 802x^2 - 800xy + 200y^2 > 0 $$
% 
% Finally, we have
%
% $$ 2x^2 + 200(2x - y)^2 > 0,$$
%
% that obviously holds for all $(x, y) \neq (0,0)$. 

%% 
% e) Analytical contours
%
% $$ f(x, y) = 100(y-x^2)^2 + (1-x)^2 $$
%
% We compute the contour f(x, y) = c, with c > 0
%
% $$ 100(y-x^2)^2 + (1-x)^2 = c $$
%
% $$ 100(y-x^2)^2 = c - (1-x)^2$$
%
% $$ 100(y-x^2) = \pm\sqrt{c - (1-x)^2}$$
% 
% $$ y = x^2 \pm \frac{\sqrt{c - (1-x)^2}}{100}$$

%% 
% f) Finite difference implementation of the gradient

% Construct central finite difference to be applied columnwise
D = 1/(2*h)*spdiags([-ones(n+1,1) ones(n+1,1)], [-1 1], n+1,n+1);
D(1,1:2) = 1/h*[-1 1];
D(n+1,n:n+1) = 1/h*[-1 1];

% Gradient via finate differences 
% d/dy
Dyf = D*f(alpha*X, alpha*Y); 
% d/dx, central finite difference applied rowwise 
Dxf = (D*(f(alpha*X, alpha*Y).')).';

% Hessian via finite difference
% d^2/dxdy
Dxyf = D*Dxf;
% d^2/dx^2
Dxxf = (D*(Dxf.')).';
% d^2/dy^2
Dyyf = D*Dyf;
% d^2/dydx
Dyxf = (D*(Dyf.')).';




