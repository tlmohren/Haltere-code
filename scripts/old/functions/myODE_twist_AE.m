function dy = myODE_twist(t,y)
 dy(1,1) =y(3);
 dy(2,1) =y(4);
 dy(3,1) =((((y(4)*sin(y(1))*sin(y(2)))/10 - 10*sin((pi*sin(400*pi*t))/2)*sin(y(1)) + sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))*(10*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + cos((pi*sin(400*pi*t))/2)*sin(y(2)) - (y(4)*cos(y(1))*sin(y(2)))/10 + sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 - ((10*sin((pi*sin(400*pi*t))/2)*sin(y(1)) + (y(4)*sin(y(1))*sin(y(2)))/10 + sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))*(10*sin((pi*sin(400*pi*t))/2)*cos(y(1)) - cos((pi*sin(400*pi*t))/2)*sin(y(2)) + (y(4)*cos(y(1))*sin(y(2)))/10 - sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 - ((cos((pi*sin(400*pi*t))/2)*sin(y(1)) - (cos((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))/10)*(y(3)^2*cos((pi*sin(400*pi*t))/2)*cos(y(1)) - (y(4)^2*sin((pi*sin(400*pi*t))/2)*sin(y(2)))/10 - 80000*pi^3*sin(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1)) - 8000*pi^3*sin(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(2)) + 40000*pi^4*cos(400*pi*t)^2*cos((pi*sin(400*pi*t))/2)*cos(y(1)) - 4000*pi^4*cos(400*pi*t)^2*sin((pi*sin(400*pi*t))/2)*sin(y(2)) + (y(3)^2*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 + (y(4)^2*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 - 8000*pi^3*sin(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)) + 4000*pi^4*cos(400*pi*t)^2*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)) + 40*pi^2*y(4)*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2)) - 400*pi^2*y(3)*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(1)) + (y(3)*y(4)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/5 + 40*pi^2*y(3)*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)) - 40*pi^2*y(4)*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2))))/10 + ((cos((pi*sin(400*pi*t))/2)*sin(y(1)) + (cos((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))/10)*(80000*pi^3*sin(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1)) - (y(4)^2*sin((pi*sin(400*pi*t))/2)*sin(y(2)))/10 - y(3)^2*cos((pi*sin(400*pi*t))/2)*cos(y(1)) - 8000*pi^3*sin(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(2)) - 40000*pi^4*cos(400*pi*t)^2*cos((pi*sin(400*pi*t))/2)*cos(y(1)) - 4000*pi^4*cos(400*pi*t)^2*sin((pi*sin(400*pi*t))/2)*sin(y(2)) + (y(3)^2*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 + (y(4)^2*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 - 8000*pi^3*sin(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)) + 4000*pi^4*cos(400*pi*t)^2*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)) + 40*pi^2*y(4)*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2)) + 400*pi^2*y(3)*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(1)) + (y(3)*y(4)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/5 + 40*pi^2*y(3)*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)) - 40*pi^2*y(4)*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2))))/10 - (3391170250170917*cos(y(1))*sin(y(1)))/2147483648 - 80*pi*y(3) - ((y(3)*cos((pi*sin(400*pi*t))/2)*cos(y(1)) - 200*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(1)) + (y(3)*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 + (y(4)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/10 + 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))*(y(3)*cos((pi*sin(400*pi*t))/2)*sin(y(1)) + (y(4)*sin((pi*sin(400*pi*t))/2)*cos(y(2)))/10 + 200*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(2)) - (y(3)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))/10 + (y(4)*cos((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))/10 + 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 - ((200*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(1)) - y(3)*cos((pi*sin(400*pi*t))/2)*cos(y(1)) + (y(3)*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 + (y(4)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/10 + 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))*((y(4)*sin((pi*sin(400*pi*t))/2)*cos(y(2)))/10 - y(3)*cos((pi*sin(400*pi*t))/2)*sin(y(1)) - 200*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(2)) - (y(3)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))/10 + (y(4)*cos((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))/10 + 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 + ((cos(y(1)) + (cos(y(2))*sin(y(1)))/10)*(y(3)^2*sin(y(1)) + y(4)*cos((pi*sin(400*pi*t))/2)*cos(y(2)) - 10*y(3)*sin((pi*sin(400*pi*t))/2)*sin(y(1)) - (y(3)^2*cos(y(1))*cos(y(2)))/10 - (y(4)^2*cos(y(1))*cos(y(2)))/10 + 2000*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1)) + (y(3)*y(4)*sin(y(1))*sin(y(2)))/5 - 200*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(2)) + y(3)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)) - y(4)*sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)) + 200*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 + ((cos(y(1)) - (cos(y(2))*sin(y(1)))/10)*(y(3)^2*sin(y(1)) - y(4)*cos((pi*sin(400*pi*t))/2)*cos(y(2)) - 10*y(3)*sin((pi*sin(400*pi*t))/2)*sin(y(1)) + (y(3)^2*cos(y(1))*cos(y(2)))/10 + (y(4)^2*cos(y(1))*cos(y(2)))/10 + 2000*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1)) - (y(3)*y(4)*sin(y(1))*sin(y(2)))/5 + 200*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(2)) - y(3)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)) + y(4)*sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)) - 200*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 - ((sin((pi*sin(400*pi*t))/2)*sin(y(1)) - (sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))/10)*(y(3)^2*sin((pi*sin(400*pi*t))/2)*cos(y(1)) - y(3)*cos(y(2))*sin(y(1)) - y(4)*cos(y(1))*sin(y(2)) - 10*y(3)*cos(y(1)) + (y(4)^2*cos((pi*sin(400*pi*t))/2)*sin(y(2)))/10 + 80000*pi^3*sin(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1)) - 8000*pi^3*sin(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(2)) + 40000*pi^4*cos(400*pi*t)^2*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + 4000*pi^4*cos(400*pi*t)^2*cos((pi*sin(400*pi*t))/2)*sin(y(2)) + (y(3)^2*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 + (y(4)^2*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 + 8000*pi^3*sin(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)) + 4000*pi^4*cos(400*pi*t)^2*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)) + 400*pi^2*y(3)*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(1)) + 40*pi^2*y(4)*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(2)) + (y(3)*y(4)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/5 - 40*pi^2*y(3)*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)) + 40*pi^2*y(4)*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2))))/10 + ((sin((pi*sin(400*pi*t))/2)*sin(y(1)) + (sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))/10)*(10*y(3)*cos(y(1)) - y(3)*cos(y(2))*sin(y(1)) - y(4)*cos(y(1))*sin(y(2)) - y(3)^2*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + (y(4)^2*cos((pi*sin(400*pi*t))/2)*sin(y(2)))/10 - 80000*pi^3*sin(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1)) - 8000*pi^3*sin(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(2)) - 40000*pi^4*cos(400*pi*t)^2*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + 4000*pi^4*cos(400*pi*t)^2*cos((pi*sin(400*pi*t))/2)*sin(y(2)) + (y(3)^2*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 + (y(4)^2*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 + 8000*pi^3*sin(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)) + 4000*pi^4*cos(400*pi*t)^2*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)) - 400*pi^2*y(3)*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(1)) + 40*pi^2*y(4)*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(2)) + (y(3)*y(4)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/5 - 40*pi^2*y(3)*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)) + 40*pi^2*y(4)*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2))))/10 + (4*cos(y(1))*(y(3)^2*sin(y(1)) - 10*y(3)*sin((pi*sin(400*pi*t))/2)*sin(y(1)) + 2000*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1))))/5 - ((cos(y(2))*sin(y(1)) - 10*cos(y(1)) + 200*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(1)) - (y(4)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/10 + 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))*(10*sin(y(1)) + cos(y(1))*cos(y(2)) - (y(4)*cos((pi*sin(400*pi*t))/2)*cos(y(2)))/10 + 200*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1)) + 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(2)) + (y(4)*sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))/10 - 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 + ((10*cos(y(1)) + cos(y(2))*sin(y(1)) - 200*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(1)) - (y(4)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/10 + 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))*(10*sin(y(1)) - cos(y(1))*cos(y(2)) + (y(4)*cos((pi*sin(400*pi*t))/2)*cos(y(2)))/10 + 200*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1)) - 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(2)) - (y(4)*sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))/10 + 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 + ((y(3)*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + 200*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(1)) - (y(3)*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 - (y(4)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/10 + 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))*(10*sin(y(1)) + cos(y(1))*cos(y(2)) - (y(4)*cos((pi*sin(400*pi*t))/2)*cos(y(2)))/10 - y(3)*sin((pi*sin(400*pi*t))/2)*sin(y(1)) + 200*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1)) + 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(2)) - (y(3)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))/10 + (y(4)*sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))/10 - 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 + ((y(3)*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + 200*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(1)) + (y(3)*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 + (y(4)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/10 - 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))*(10*sin(y(1)) - cos(y(1))*cos(y(2)) + (y(4)*cos((pi*sin(400*pi*t))/2)*cos(y(2)))/10 - y(3)*sin((pi*sin(400*pi*t))/2)*sin(y(1)) + 200*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1)) - 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(2)) + (y(3)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))/10 - (y(4)*sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))/10 + 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 + ((y(3)*sin(y(1)) + (y(3)*cos(y(1))*cos(y(2)))/10 - (y(4)*sin(y(1))*sin(y(2)))/10)*(y(3)*cos(y(1)) - 10*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + cos((pi*sin(400*pi*t))/2)*sin(y(2)) - (y(3)*cos(y(2))*sin(y(1)))/10 - (y(4)*cos(y(1))*sin(y(2)))/10 + sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 + ((y(3)*sin(y(1)) - (y(3)*cos(y(1))*cos(y(2)))/10 + (y(4)*sin(y(1))*sin(y(2)))/10)*(y(3)*cos(y(1)) - 10*sin((pi*sin(400*pi*t))/2)*cos(y(1)) - cos((pi*sin(400*pi*t))/2)*sin(y(2)) + (y(3)*cos(y(2))*sin(y(1)))/10 + (y(4)*cos(y(1))*sin(y(2)))/10 - sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 + (((y(4)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/10 - 200*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(1)) + 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))*((y(4)*sin((pi*sin(400*pi*t))/2)*cos(y(2)))/10 + 200*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(2)) + (y(4)*cos((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))/10 + 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 + ((200*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(1)) + (y(4)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/10 + 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))*((y(4)*sin((pi*sin(400*pi*t))/2)*cos(y(2)))/10 - 200*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(2)) + (y(4)*cos((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))/10 + 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 + (4*(10*cos(y(1)) - 200*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(1)))*(10*sin(y(1)) + 200*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1))))/5 + (3391170250170917*cos((pi*sin(400*pi*t))/2)^2*cos(y(1))*sin(y(1)))/2147483648 - 80*sin((pi*sin(400*pi*t))/2)^2*cos(y(1))*sin(y(1)) - (4*sin((pi*sin(400*pi*t))/2)*sin(y(1))*(y(3)^2*sin((pi*sin(400*pi*t))/2)*cos(y(1)) - 10*y(3)*cos(y(1)) + 80000*pi^3*sin(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1)) + 40000*pi^4*cos(400*pi*t)^2*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + 400*pi^2*y(3)*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(1))))/5 - (4*cos((pi*sin(400*pi*t))/2)*sin(y(1))*(y(3)^2*cos((pi*sin(400*pi*t))/2)*cos(y(1)) - 80000*pi^3*sin(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + 40000*pi^4*cos(400*pi*t)^2*cos((pi*sin(400*pi*t))/2)*cos(y(1)) - 400*pi^2*y(3)*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(1))))/5 + (4*y(3)*sin(y(1))*(y(3)*cos(y(1)) - 10*sin((pi*sin(400*pi*t))/2)*cos(y(1))))/5 - (4*y(3)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*(y(3)*cos((pi*sin(400*pi*t))/2)*sin(y(1)) + 200*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1))))/5 + (4*y(3)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*(10*sin(y(1)) - y(3)*sin((pi*sin(400*pi*t))/2)*sin(y(1)) + 200*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1))))/5 + 160*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(1))*(y(3)*cos((pi*sin(400*pi*t))/2)*sin(y(1)) + 200*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1))) + 160*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(1))*(10*sin(y(1)) - y(3)*sin((pi*sin(400*pi*t))/2)*sin(y(1)) + 200*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1))) - 32000*pi^4*cos(400*pi*t)^2*sin((pi*sin(400*pi*t))/2)^2*cos(y(1))*sin(y(1)))/((4*cos(y(1))^2)/5 + (cos((pi*sin(400*pi*t))/2)*sin(y(1)) - (cos((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))/10)^2/10 + (cos((pi*sin(400*pi*t))/2)*sin(y(1)) + (cos((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))/10)^2/10 + (sin((pi*sin(400*pi*t))/2)*sin(y(1)) - (sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))/10)^2/10 + (sin((pi*sin(400*pi*t))/2)*sin(y(1)) + (sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))/10)^2/10 + (4*cos((pi*sin(400*pi*t))/2)^2*sin(y(1))^2)/5 + (cos(y(1)) - (cos(y(2))*sin(y(1)))/10)^2/10 + (cos(y(1)) + (cos(y(2))*sin(y(1)))/10)^2/10 + (4*sin((pi*sin(400*pi*t))/2)^2*sin(y(1))^2)/5);
 dy(4,1) =-((3391170250170917*y(2))/4294967296 + 80*pi*y(4) + (((sin((pi*sin(400*pi*t))/2)*cos(y(2)))/10 + (cos((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))/10)*(y(3)^2*cos((pi*sin(400*pi*t))/2)*cos(y(1)) - (y(4)^2*sin((pi*sin(400*pi*t))/2)*sin(y(2)))/10 - 80000*pi^3*sin(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1)) - 8000*pi^3*sin(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(2)) + 40000*pi^4*cos(400*pi*t)^2*cos((pi*sin(400*pi*t))/2)*cos(y(1)) - 4000*pi^4*cos(400*pi*t)^2*sin((pi*sin(400*pi*t))/2)*sin(y(2)) + (y(3)^2*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 + (y(4)^2*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 - 8000*pi^3*sin(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)) + 4000*pi^4*cos(400*pi*t)^2*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)) + 40*pi^2*y(4)*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2)) - 400*pi^2*y(3)*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(1)) + (y(3)*y(4)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/5 + 40*pi^2*y(3)*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)) - 40*pi^2*y(4)*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2))))/10 + (((sin((pi*sin(400*pi*t))/2)*cos(y(2)))/10 + (cos((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))/10)*(80000*pi^3*sin(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1)) - (y(4)^2*sin((pi*sin(400*pi*t))/2)*sin(y(2)))/10 - y(3)^2*cos((pi*sin(400*pi*t))/2)*cos(y(1)) - 8000*pi^3*sin(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(2)) - 40000*pi^4*cos(400*pi*t)^2*cos((pi*sin(400*pi*t))/2)*cos(y(1)) - 4000*pi^4*cos(400*pi*t)^2*sin((pi*sin(400*pi*t))/2)*sin(y(2)) + (y(3)^2*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 + (y(4)^2*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 - 8000*pi^3*sin(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)) + 4000*pi^4*cos(400*pi*t)^2*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)) + 40*pi^2*y(4)*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2)) + 400*pi^2*y(3)*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(1)) + (y(3)*y(4)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/5 + 40*pi^2*y(3)*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)) - 40*pi^2*y(4)*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2))))/10 - (((cos((pi*sin(400*pi*t))/2)*cos(y(2)))/10 - (sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))/10)*(y(3)^2*sin((pi*sin(400*pi*t))/2)*cos(y(1)) - y(3)*cos(y(2))*sin(y(1)) - y(4)*cos(y(1))*sin(y(2)) - 10*y(3)*cos(y(1)) + (y(4)^2*cos((pi*sin(400*pi*t))/2)*sin(y(2)))/10 + 80000*pi^3*sin(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1)) - 8000*pi^3*sin(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(2)) + 40000*pi^4*cos(400*pi*t)^2*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + 4000*pi^4*cos(400*pi*t)^2*cos((pi*sin(400*pi*t))/2)*sin(y(2)) + (y(3)^2*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 + (y(4)^2*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 + 8000*pi^3*sin(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)) + 4000*pi^4*cos(400*pi*t)^2*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)) + 400*pi^2*y(3)*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(1)) + 40*pi^2*y(4)*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(2)) + (y(3)*y(4)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/5 - 40*pi^2*y(3)*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)) + 40*pi^2*y(4)*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2))))/10 - (((cos((pi*sin(400*pi*t))/2)*cos(y(2)))/10 - (sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))/10)*(10*y(3)*cos(y(1)) - y(3)*cos(y(2))*sin(y(1)) - y(4)*cos(y(1))*sin(y(2)) - y(3)^2*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + (y(4)^2*cos((pi*sin(400*pi*t))/2)*sin(y(2)))/10 - 80000*pi^3*sin(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1)) - 8000*pi^3*sin(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(2)) - 40000*pi^4*cos(400*pi*t)^2*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + 4000*pi^4*cos(400*pi*t)^2*cos((pi*sin(400*pi*t))/2)*sin(y(2)) + (y(3)^2*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 + (y(4)^2*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 + 8000*pi^3*sin(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)) + 4000*pi^4*cos(400*pi*t)^2*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)) - 400*pi^2*y(3)*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(1)) + 40*pi^2*y(4)*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(2)) + (y(3)*y(4)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/5 - 40*pi^2*y(3)*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)) + 40*pi^2*y(4)*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2))))/10 + ((20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2)) - (y(4)*sin((pi*sin(400*pi*t))/2)*sin(y(2)))/10 + (y(3)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/10 + (y(4)*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 - 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))*(y(3)*cos((pi*sin(400*pi*t))/2)*sin(y(1)) + (y(4)*sin((pi*sin(400*pi*t))/2)*cos(y(2)))/10 + 200*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(2)) - (y(3)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))/10 + (y(4)*cos((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))/10 + 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 + ((20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2)) - (y(4)*sin((pi*sin(400*pi*t))/2)*sin(y(2)))/10 + (y(3)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/10 + (y(4)*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 - 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))*((y(4)*sin((pi*sin(400*pi*t))/2)*cos(y(2)))/10 - y(3)*cos((pi*sin(400*pi*t))/2)*sin(y(1)) - 200*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(2)) - (y(3)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))/10 + (y(4)*cos((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))/10 + 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 - ((20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(2)) - cos(y(1))*sin(y(2)) + (y(3)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/10 + 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))*(10*sin(y(1)) + cos(y(1))*cos(y(2)) - y(3)*sin((pi*sin(400*pi*t))/2)*sin(y(1)) + 200*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1)) + 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(2)) - (y(3)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))/10 - 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 + ((20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(2)) - cos(y(1))*sin(y(2)) + (y(3)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/10 + 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))*(10*sin(y(1)) - cos(y(1))*cos(y(2)) - y(3)*sin((pi*sin(400*pi*t))/2)*sin(y(1)) + 200*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1)) - 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(2)) + (y(3)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))/10 + 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 - ((cos((pi*sin(400*pi*t))/2)*cos(y(2)) + (y(3)*sin(y(1))*sin(y(2)))/10 - sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))*(y(3)*cos(y(1)) - 10*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + cos((pi*sin(400*pi*t))/2)*sin(y(2)) - (y(3)*cos(y(2))*sin(y(1)))/10 + sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 - ((cos((pi*sin(400*pi*t))/2)*cos(y(2)) + (y(3)*sin(y(1))*sin(y(2)))/10 - sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))*(10*sin((pi*sin(400*pi*t))/2)*cos(y(1)) - y(3)*cos(y(1)) + cos((pi*sin(400*pi*t))/2)*sin(y(2)) - (y(3)*cos(y(2))*sin(y(1)))/10 + sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 + ((20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2)) + (y(3)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/10 - 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))*(y(3)*cos((pi*sin(400*pi*t))/2)*sin(y(1)) + 200*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1)) - 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(2)) + (y(3)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))/10 - 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 - ((20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2)) + (y(3)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/10 - 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))*(y(3)*cos((pi*sin(400*pi*t))/2)*sin(y(1)) + 200*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(2)) - (y(3)*cos((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))/10 + 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 + (((y(4)*cos((pi*sin(400*pi*t))/2)*sin(y(2)))/10 + 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(2)) + (y(3)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/10 + (y(4)*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 + 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))*(10*sin(y(1)) + cos(y(1))*cos(y(2)) - (y(4)*cos((pi*sin(400*pi*t))/2)*cos(y(2)))/10 - y(3)*sin((pi*sin(400*pi*t))/2)*sin(y(1)) + 200*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1)) + 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(2)) - (y(3)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))/10 + (y(4)*sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))/10 - 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 - (((y(4)*cos((pi*sin(400*pi*t))/2)*sin(y(2)))/10 + 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*cos(y(2)) + (y(3)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*sin(y(2)))/10 + (y(4)*sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1)))/10 + 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))*(10*sin(y(1)) - cos(y(1))*cos(y(2)) + (y(4)*cos((pi*sin(400*pi*t))/2)*cos(y(2)))/10 - y(3)*sin((pi*sin(400*pi*t))/2)*sin(y(1)) + 200*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1)) - 20*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(2)) + (y(3)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)))/10 - (y(4)*sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))/10 + 20*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/10 - (cos(y(1))*sin(y(2))*(y(3)^2*sin(y(1)) + y(4)*cos((pi*sin(400*pi*t))/2)*cos(y(2)) - 10*y(3)*sin((pi*sin(400*pi*t))/2)*sin(y(1)) - (y(3)^2*cos(y(1))*cos(y(2)))/10 - (y(4)^2*cos(y(1))*cos(y(2)))/10 + 2000*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1)) + (y(3)*y(4)*sin(y(1))*sin(y(2)))/5 - 200*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(2)) + y(3)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)) - y(4)*sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)) + 200*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/100 + (cos(y(1))*sin(y(2))*(y(3)^2*sin(y(1)) - y(4)*cos((pi*sin(400*pi*t))/2)*cos(y(2)) - 10*y(3)*sin((pi*sin(400*pi*t))/2)*sin(y(1)) + (y(3)^2*cos(y(1))*cos(y(2)))/10 + (y(4)^2*cos(y(1))*cos(y(2)))/10 + 2000*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(1)) - (y(3)*y(4)*sin(y(1))*sin(y(2)))/5 + 200*pi^2*cos(400*pi*t)*sin((pi*sin(400*pi*t))/2)*sin(y(2)) - y(3)*sin((pi*sin(400*pi*t))/2)*cos(y(1))*cos(y(2)) + y(4)*sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)) - 200*pi^2*cos(400*pi*t)*cos((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/100 - (y(4)*cos(y(1))*cos(y(2))*(y(3)*cos(y(1)) - 10*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + cos((pi*sin(400*pi*t))/2)*sin(y(2)) - (y(3)*cos(y(2))*sin(y(1)))/10 - (y(4)*cos(y(1))*sin(y(2)))/10 + sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/100 + (y(4)*cos(y(1))*cos(y(2))*(y(3)*cos(y(1)) - 10*sin((pi*sin(400*pi*t))/2)*cos(y(1)) - cos((pi*sin(400*pi*t))/2)*sin(y(2)) + (y(3)*cos(y(2))*sin(y(1)))/10 + (y(4)*cos(y(1))*sin(y(2)))/10 - sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/100 + (y(3)*sin(y(1))*sin(y(2))*(y(3)*cos(y(1)) - 10*sin((pi*sin(400*pi*t))/2)*cos(y(1)) + cos((pi*sin(400*pi*t))/2)*sin(y(2)) - (y(3)*cos(y(2))*sin(y(1)))/10 - (y(4)*cos(y(1))*sin(y(2)))/10 + sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/100 - (y(3)*sin(y(1))*sin(y(2))*(y(3)*cos(y(1)) - 10*sin((pi*sin(400*pi*t))/2)*cos(y(1)) - cos((pi*sin(400*pi*t))/2)*sin(y(2)) + (y(3)*cos(y(2))*sin(y(1)))/10 + (y(4)*cos(y(1))*sin(y(2)))/10 - sin((pi*sin(400*pi*t))/2)*cos(y(2))*sin(y(1))))/100)/(((cos((pi*sin(400*pi*t))/2)*cos(y(2)))/10 - (sin((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))/10)^2/5 + ((sin((pi*sin(400*pi*t))/2)*cos(y(2)))/10 + (cos((pi*sin(400*pi*t))/2)*sin(y(1))*sin(y(2)))/10)^2/5 + (cos(y(1))^2*sin(y(2))^2)/500);
 t