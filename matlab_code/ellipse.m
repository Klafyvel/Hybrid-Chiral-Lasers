% Compute the ellipse associated to a polarization state

function [X,Y] = ellipse(epsilon, theta, x0, y0, a, n)
    if nargin < 6
        n = 200;
    end
    if nargin < 5
        a = 1;
    end
    if nargin < 3
        x0 = 0;
        y0 = 0;
    end
    b = a*tan(epsilon);
    t = linspace(0,2*pi,n);
    X = a*cos(t);
    Y = b*sin(t);
    v = rotation(theta) * [X;Y];
    X = v(1,:) + x0;
    Y = v(2,:) + y0;
end
