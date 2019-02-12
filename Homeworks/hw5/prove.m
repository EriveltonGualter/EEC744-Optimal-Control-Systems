clear
syms xf yf thetaf theta 

x = xf + xf + yf/(2*cos(thetaf)^2) * (2*(thetaf-theta) + sin(2*thetaf) - sin(2*theta));
y = yf*cos(theta)^2/(cos(thetaf)^2);

ydot = diff(y,theta)*1/(diff(x,theta));

yddot = diff(ydot,theta)*1/(diff(x,theta));

simplify(ydot)

simplify(y*yddot)

simplify(1 + ydot^2 + 2*y*yddot)


