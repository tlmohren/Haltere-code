function dy = myODE_twist(t,y)
 dy(1,1) =y(3);
 dy(2,1) =y(4);
 dy(3,1) =-(((cos((pi*sin(400*t*pi))/2)*sin(y(1)) - (cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/100)*(10*y(3)*cos(y(1)) + (y(3)*cos(y(2))*sin(y(1)))/10 + (y(4)*cos(y(1))*sin(y(2)))/10 + y(3)^2*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - (y(4)^2*sin((pi*sin(400*t*pi))/2)*sin(y(2)))/100 - 80000*pi^3*sin(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1)) - 800*pi^3*sin(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(2)) + 40000*pi^4*cos(400*t*pi)^2*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - 400*pi^4*cos(400*t*pi)^2*sin((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(3)^2*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)^2*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 - 800*pi^3*sin(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)) + 400*pi^4*cos(400*t*pi)^2*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)) + 4*pi^2*y(4)*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(2)) - 400*pi^2*y(3)*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(1)) + (y(3)*y(4)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/50 + 4*pi^2*y(3)*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)) - 4*pi^2*y(4)*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2))))/20 - (9*(10*cos(y(1)) - 200*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(1)))*(10*sin(y(1)) + 200*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1))))/10 - ((cos((pi*sin(400*t*pi))/2)*sin(y(1)) + (cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/100)*((y(3)*cos(y(2))*sin(y(1)))/10 - 10*y(3)*cos(y(1)) + (y(4)*cos(y(1))*sin(y(2)))/10 - y(3)^2*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - (y(4)^2*sin((pi*sin(400*t*pi))/2)*sin(y(2)))/100 + 80000*pi^3*sin(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1)) - 800*pi^3*sin(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(2)) - 40000*pi^4*cos(400*t*pi)^2*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - 400*pi^4*cos(400*t*pi)^2*sin((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(3)^2*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)^2*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 - 800*pi^3*sin(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)) + 400*pi^4*cos(400*t*pi)^2*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)) + 4*pi^2*y(4)*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(2)) + 400*pi^2*y(3)*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(1)) + (y(3)*y(4)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/50 + 4*pi^2*y(3)*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)) - 4*pi^2*y(4)*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2))))/20 + (2649351757946029*cos(y(1))*sin(y(1)))/16777216 + 80*pi*y(3) - ((cos(y(1)) - (cos(y(2))*sin(y(1)))/100)*(y(3)^2*sin(y(1)) + 10*y(3)*cos((pi*sin(400*t*pi))/2)*sin(y(1)) - (y(4)*sin((pi*sin(400*t*pi))/2)*cos(y(2)))/10 + (y(3)^2*cos(y(1))*cos(y(2)))/100 + (y(4)^2*cos(y(1))*cos(y(2)))/100 - (y(3)*y(4)*sin(y(1))*sin(y(2)))/50 + 2000*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1)) - 20*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(3)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/10 - (y(4)*cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/10 - 20*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1))))/20 - ((cos(y(1)) + (cos(y(2))*sin(y(1)))/100)*(y(3)^2*sin(y(1)) + 10*y(3)*cos((pi*sin(400*t*pi))/2)*sin(y(1)) + (y(4)*sin((pi*sin(400*t*pi))/2)*cos(y(2)))/10 - (y(3)^2*cos(y(1))*cos(y(2)))/100 - (y(4)^2*cos(y(1))*cos(y(2)))/100 + (y(3)*y(4)*sin(y(1))*sin(y(2)))/50 + 2000*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1)) + 20*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(2)) - (y(3)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/10 + (y(4)*cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/10 + 20*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1))))/20 + ((sin((pi*sin(400*t*pi))/2)*sin(y(1)) - (sin((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/100)*(y(3)^2*sin((pi*sin(400*t*pi))/2)*cos(y(1)) + (y(4)^2*cos((pi*sin(400*t*pi))/2)*sin(y(2)))/100 + 80000*pi^3*sin(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - 800*pi^3*sin(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(2)) + 40000*pi^4*cos(400*t*pi)^2*sin((pi*sin(400*t*pi))/2)*cos(y(1)) + 400*pi^4*cos(400*t*pi)^2*cos((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(3)^2*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)^2*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + 800*pi^3*sin(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)) + 400*pi^4*cos(400*t*pi)^2*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)) + 400*pi^2*y(3)*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(1)) + 4*pi^2*y(4)*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2)) + (y(3)*y(4)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/50 - 4*pi^2*y(3)*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)) + 4*pi^2*y(4)*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2))))/20 - ((sin((pi*sin(400*t*pi))/2)*sin(y(1)) + (sin((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/100)*((y(4)^2*cos((pi*sin(400*t*pi))/2)*sin(y(2)))/100 - y(3)^2*sin((pi*sin(400*t*pi))/2)*cos(y(1)) - 80000*pi^3*sin(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - 800*pi^3*sin(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(2)) - 40000*pi^4*cos(400*t*pi)^2*sin((pi*sin(400*t*pi))/2)*cos(y(1)) + 400*pi^4*cos(400*t*pi)^2*cos((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(3)^2*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)^2*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + 800*pi^3*sin(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)) + 400*pi^4*cos(400*t*pi)^2*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)) - 400*pi^2*y(3)*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(1)) + 4*pi^2*y(4)*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2)) + (y(3)*y(4)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/50 - 4*pi^2*y(3)*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)) + 4*pi^2*y(4)*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2))))/20 - (9*cos(y(1))*(y(3)^2*sin(y(1)) + 10*y(3)*cos((pi*sin(400*t*pi))/2)*sin(y(1)) + 2000*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1))))/10 + ((y(3)*sin((pi*sin(400*t*pi))/2)*cos(y(1)) + 200*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(1)) - (y(3)*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 - (y(4)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))*((y(4)*cos((pi*sin(400*t*pi))/2)*cos(y(2)))/100 + y(3)*sin((pi*sin(400*t*pi))/2)*sin(y(1)) - 200*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(3)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/100 - (y(4)*sin((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/100 + 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1))))/20 - ((y(3)*sin((pi*sin(400*t*pi))/2)*cos(y(1)) + 200*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(1)) + (y(3)*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/100 - 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))*((y(4)*cos((pi*sin(400*t*pi))/2)*cos(y(2)))/100 - y(3)*sin((pi*sin(400*t*pi))/2)*sin(y(1)) + 200*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(3)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/100 - (y(4)*sin((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/100 + 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1))))/20 + ((y(3)*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - 200*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(1)) + (y(3)*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))*(10*sin(y(1)) - (cos(y(1))*cos(y(2)))/10 + y(3)*cos((pi*sin(400*t*pi))/2)*sin(y(1)) + (y(4)*sin((pi*sin(400*t*pi))/2)*cos(y(2)))/100 + 200*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1)) + 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(2)) - (y(3)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/100 + (y(4)*cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/100 + 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1))))/20 - ((200*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(1)) - y(3)*cos((pi*sin(400*t*pi))/2)*cos(y(1)) + (y(3)*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))*(10*sin(y(1)) + (cos(y(1))*cos(y(2)))/10 + y(3)*cos((pi*sin(400*t*pi))/2)*sin(y(1)) - (y(4)*sin((pi*sin(400*t*pi))/2)*cos(y(2)))/100 + 200*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1)) - 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(3)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/100 - (y(4)*cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/100 - 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1))))/20 - (((cos(y(2))*sin(y(1)))/10 - 10*cos(y(1)) + 200*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(1)) + (y(4)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))*((y(4)*sin((pi*sin(400*t*pi))/2)*cos(y(2)))/100 - (cos(y(1))*cos(y(2)))/10 - 10*sin(y(1)) - 200*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1)) + 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(4)*cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/100 + 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1))))/20 - ((10*cos(y(1)) + (cos(y(2))*sin(y(1)))/10 - 200*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(1)) + (y(4)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))*(10*sin(y(1)) - (cos(y(1))*cos(y(2)))/10 + (y(4)*sin((pi*sin(400*t*pi))/2)*cos(y(2)))/100 + 200*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1)) + 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(4)*cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/100 + 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1))))/20 + ((200*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(1)) - (y(4)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))*(200*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - (y(4)*cos((pi*sin(400*t*pi))/2)*cos(y(2)))/100 + 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(4)*sin((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/100 - 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1))))/20 + ((200*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(1)) + (y(4)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/100 - 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))*((y(4)*cos((pi*sin(400*t*pi))/2)*cos(y(2)))/100 + 200*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(2)) - (y(4)*sin((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/100 + 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1))))/20 - ((y(3)*sin(y(1)) - (y(3)*cos(y(1))*cos(y(2)))/100 + (y(4)*sin(y(1))*sin(y(2)))/100)*(y(3)*cos(y(1)) - (sin((pi*sin(400*t*pi))/2)*sin(y(2)))/10 + 10*cos((pi*sin(400*t*pi))/2)*cos(y(1)) + (y(3)*cos(y(2))*sin(y(1)))/100 + (y(4)*cos(y(1))*sin(y(2)))/100 + (cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/10))/20 - ((y(3)*sin(y(1)) + (y(3)*cos(y(1))*cos(y(2)))/100 - (y(4)*sin(y(1))*sin(y(2)))/100)*((sin((pi*sin(400*t*pi))/2)*sin(y(2)))/10 + y(3)*cos(y(1)) + 10*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - (y(3)*cos(y(2))*sin(y(1)))/100 - (y(4)*cos(y(1))*sin(y(2)))/100 - (cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/10))/20 + ((10*cos((pi*sin(400*t*pi))/2)*sin(y(1)) - (y(4)*sin(y(1))*sin(y(2)))/100 + (cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/10)*((sin((pi*sin(400*t*pi))/2)*sin(y(2)))/10 + 10*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - (y(4)*cos(y(1))*sin(y(2)))/100 - (cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/10))/20 + ((10*cos((pi*sin(400*t*pi))/2)*sin(y(1)) + (y(4)*sin(y(1))*sin(y(2)))/100 - (cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/10)*(10*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - (sin((pi*sin(400*t*pi))/2)*sin(y(2)))/10 + (y(4)*cos(y(1))*sin(y(2)))/100 + (cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/10))/20 - (2649350247996589*cos((pi*sin(400*t*pi))/2)^2*cos(y(1))*sin(y(1)))/16777216 + (9*cos((pi*sin(400*t*pi))/2)*sin(y(1))*(10*y(3)*cos(y(1)) + y(3)^2*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - 80000*pi^3*sin(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1)) + 40000*pi^4*cos(400*t*pi)^2*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - 400*pi^2*y(3)*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(1))))/10 + (9*sin((pi*sin(400*t*pi))/2)*sin(y(1))*(y(3)^2*sin((pi*sin(400*t*pi))/2)*cos(y(1)) + 80000*pi^3*sin(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1)) + 40000*pi^4*cos(400*t*pi)^2*sin((pi*sin(400*t*pi))/2)*cos(y(1)) + 400*pi^2*y(3)*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(1))))/10 - (9*y(3)*sin(y(1))*(y(3)*cos(y(1)) + 10*cos((pi*sin(400*t*pi))/2)*cos(y(1))))/10 + (9*y(3)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*(y(3)*sin((pi*sin(400*t*pi))/2)*sin(y(1)) - 200*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1))))/10 + (9*y(3)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*(10*sin(y(1)) + y(3)*cos((pi*sin(400*t*pi))/2)*sin(y(1)) + 200*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1))))/10 + 180*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(1))*(y(3)*sin((pi*sin(400*t*pi))/2)*sin(y(1)) - 200*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1))) - 180*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(1))*(10*sin(y(1)) + y(3)*cos((pi*sin(400*t*pi))/2)*sin(y(1)) + 200*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1))) + 36000*pi^4*cos(400*t*pi)^2*cos((pi*sin(400*t*pi))/2)^2*cos(y(1))*sin(y(1)))/((9*cos(y(1))^2)/10 + (cos((pi*sin(400*t*pi))/2)*sin(y(1)) - (cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/100)^2/20 + (cos((pi*sin(400*t*pi))/2)*sin(y(1)) + (cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/100)^2/20 + (sin((pi*sin(400*t*pi))/2)*sin(y(1)) - (sin((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/100)^2/20 + (sin((pi*sin(400*t*pi))/2)*sin(y(1)) + (sin((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/100)^2/20 + (9*cos((pi*sin(400*t*pi))/2)^2*sin(y(1))^2)/10 + (cos(y(1)) - (cos(y(2))*sin(y(1)))/100)^2/20 + (cos(y(1)) + (cos(y(2))*sin(y(1)))/100)^2/20 + (9*sin((pi*sin(400*t*pi))/2)^2*sin(y(1))^2)/10);
 dy(4,1) =-((2649351757946029*y(2))/33554432 + 80*pi*y(4) + (((sin((pi*sin(400*t*pi))/2)*cos(y(2)))/100 + (cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/100)*(10*y(3)*cos(y(1)) + (y(3)*cos(y(2))*sin(y(1)))/10 + (y(4)*cos(y(1))*sin(y(2)))/10 + y(3)^2*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - (y(4)^2*sin((pi*sin(400*t*pi))/2)*sin(y(2)))/100 - 80000*pi^3*sin(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1)) - 800*pi^3*sin(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(2)) + 40000*pi^4*cos(400*t*pi)^2*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - 400*pi^4*cos(400*t*pi)^2*sin((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(3)^2*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)^2*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 - 800*pi^3*sin(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)) + 400*pi^4*cos(400*t*pi)^2*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)) + 4*pi^2*y(4)*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(2)) - 400*pi^2*y(3)*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(1)) + (y(3)*y(4)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/50 + 4*pi^2*y(3)*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)) - 4*pi^2*y(4)*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2))))/20 + (((sin((pi*sin(400*t*pi))/2)*cos(y(2)))/100 + (cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/100)*((y(3)*cos(y(2))*sin(y(1)))/10 - 10*y(3)*cos(y(1)) + (y(4)*cos(y(1))*sin(y(2)))/10 - y(3)^2*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - (y(4)^2*sin((pi*sin(400*t*pi))/2)*sin(y(2)))/100 + 80000*pi^3*sin(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1)) - 800*pi^3*sin(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(2)) - 40000*pi^4*cos(400*t*pi)^2*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - 400*pi^4*cos(400*t*pi)^2*sin((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(3)^2*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)^2*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 - 800*pi^3*sin(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)) + 400*pi^4*cos(400*t*pi)^2*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)) + 4*pi^2*y(4)*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(2)) + 400*pi^2*y(3)*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(1)) + (y(3)*y(4)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/50 + 4*pi^2*y(3)*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)) - 4*pi^2*y(4)*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2))))/20 - (((cos((pi*sin(400*t*pi))/2)*cos(y(2)))/100 - (sin((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/100)*(y(3)^2*sin((pi*sin(400*t*pi))/2)*cos(y(1)) + (y(4)^2*cos((pi*sin(400*t*pi))/2)*sin(y(2)))/100 + 80000*pi^3*sin(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - 800*pi^3*sin(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(2)) + 40000*pi^4*cos(400*t*pi)^2*sin((pi*sin(400*t*pi))/2)*cos(y(1)) + 400*pi^4*cos(400*t*pi)^2*cos((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(3)^2*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)^2*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + 800*pi^3*sin(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)) + 400*pi^4*cos(400*t*pi)^2*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)) + 400*pi^2*y(3)*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(1)) + 4*pi^2*y(4)*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2)) + (y(3)*y(4)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/50 - 4*pi^2*y(3)*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)) + 4*pi^2*y(4)*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2))))/20 - (((cos((pi*sin(400*t*pi))/2)*cos(y(2)))/100 - (sin((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/100)*((y(4)^2*cos((pi*sin(400*t*pi))/2)*sin(y(2)))/100 - y(3)^2*sin((pi*sin(400*t*pi))/2)*cos(y(1)) - 80000*pi^3*sin(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - 800*pi^3*sin(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(2)) - 40000*pi^4*cos(400*t*pi)^2*sin((pi*sin(400*t*pi))/2)*cos(y(1)) + 400*pi^4*cos(400*t*pi)^2*cos((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(3)^2*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)^2*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + 800*pi^3*sin(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)) + 400*pi^4*cos(400*t*pi)^2*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)) - 400*pi^2*y(3)*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(1)) + 4*pi^2*y(4)*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2)) + (y(3)*y(4)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/50 - 4*pi^2*y(3)*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)) + 4*pi^2*y(4)*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2))))/20 - (((cos(y(1))*sin(y(2)))/10 + 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(2)) + (y(3)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/100 - 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))*(10*sin(y(1)) - (cos(y(1))*cos(y(2)))/10 + y(3)*cos((pi*sin(400*t*pi))/2)*sin(y(1)) + 200*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1)) + 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(2)) - (y(3)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/100 + 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1))))/20 + (((cos(y(1))*sin(y(2)))/10 + 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(2)) + (y(3)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/100 - 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))*(10*sin(y(1)) + (cos(y(1))*cos(y(2)))/10 + y(3)*cos((pi*sin(400*t*pi))/2)*sin(y(1)) + 200*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1)) - 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(3)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/100 - 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1))))/20 - (((y(4)*cos((pi*sin(400*t*pi))/2)*sin(y(2)))/100 + 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2)) + (y(3)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + (y(4)*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))*((y(4)*cos((pi*sin(400*t*pi))/2)*cos(y(2)))/100 + y(3)*sin((pi*sin(400*t*pi))/2)*sin(y(1)) - 200*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(3)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/100 - (y(4)*sin((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/100 + 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1))))/20 - (((y(4)*cos((pi*sin(400*t*pi))/2)*sin(y(2)))/100 + 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2)) + (y(3)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + (y(4)*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))*((y(4)*cos((pi*sin(400*t*pi))/2)*cos(y(2)))/100 - y(3)*sin((pi*sin(400*t*pi))/2)*sin(y(1)) + 200*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(3)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/100 - (y(4)*sin((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/100 + 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1))))/20 + (((sin((pi*sin(400*t*pi))/2)*cos(y(2)))/10 + (y(3)*sin(y(1))*sin(y(2)))/100 + (cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/10)*(y(3)*cos(y(1)) - (sin((pi*sin(400*t*pi))/2)*sin(y(2)))/10 + 10*cos((pi*sin(400*t*pi))/2)*cos(y(1)) + (y(3)*cos(y(2))*sin(y(1)))/100 + (cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/10))/20 - (((sin((pi*sin(400*t*pi))/2)*cos(y(2)))/10 + (y(3)*sin(y(1))*sin(y(2)))/100 + (cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/10)*((sin((pi*sin(400*t*pi))/2)*sin(y(2)))/10 + y(3)*cos(y(1)) + 10*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - (y(3)*cos(y(2))*sin(y(1)))/100 - (cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/10))/20 + ((2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(2)) - (y(4)*sin((pi*sin(400*t*pi))/2)*sin(y(2)))/100 + (y(3)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + (y(4)*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 - 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))*(10*sin(y(1)) - (cos(y(1))*cos(y(2)))/10 + y(3)*cos((pi*sin(400*t*pi))/2)*sin(y(1)) + (y(4)*sin((pi*sin(400*t*pi))/2)*cos(y(2)))/100 + 200*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1)) + 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(2)) - (y(3)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/100 + (y(4)*cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/100 + 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1))))/20 - ((2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(2)) - (y(4)*sin((pi*sin(400*t*pi))/2)*sin(y(2)))/100 + (y(3)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + (y(4)*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/100 - 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))*(10*sin(y(1)) + (cos(y(1))*cos(y(2)))/10 + y(3)*cos((pi*sin(400*t*pi))/2)*sin(y(1)) - (y(4)*sin((pi*sin(400*t*pi))/2)*cos(y(2)))/100 + 200*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1)) - 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(3)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/100 - (y(4)*cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/100 - 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1))))/20 + ((2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2)) + (y(3)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))*(y(3)*sin((pi*sin(400*t*pi))/2)*sin(y(1)) - 200*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(3)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/100 + 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1))))/20 + ((2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2)) + (y(3)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))*(200*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - y(3)*sin((pi*sin(400*t*pi))/2)*sin(y(1)) - 2*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(3)*sin((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/100 + 2*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1))))/20 + (cos(y(1))*sin(y(2))*(y(3)^2*sin(y(1)) + 10*y(3)*cos((pi*sin(400*t*pi))/2)*sin(y(1)) - (y(4)*sin((pi*sin(400*t*pi))/2)*cos(y(2)))/10 + (y(3)^2*cos(y(1))*cos(y(2)))/100 + (y(4)^2*cos(y(1))*cos(y(2)))/100 - (y(3)*y(4)*sin(y(1))*sin(y(2)))/50 + 2000*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1)) - 20*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(2)) + (y(3)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/10 - (y(4)*cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/10 - 20*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1))))/2000 - (cos(y(1))*sin(y(2))*(y(3)^2*sin(y(1)) + 10*y(3)*cos((pi*sin(400*t*pi))/2)*sin(y(1)) + (y(4)*sin((pi*sin(400*t*pi))/2)*cos(y(2)))/10 - (y(3)^2*cos(y(1))*cos(y(2)))/100 - (y(4)^2*cos(y(1))*cos(y(2)))/100 + (y(3)*y(4)*sin(y(1))*sin(y(2)))/50 + 2000*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(1)) + 20*pi^2*cos(400*t*pi)*cos((pi*sin(400*t*pi))/2)*sin(y(2)) - (y(3)*cos((pi*sin(400*t*pi))/2)*cos(y(1))*cos(y(2)))/10 + (y(4)*cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/10 + 20*pi^2*cos(400*t*pi)*sin((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1))))/2000 + (y(4)*cos(y(1))*cos(y(2))*(y(3)*cos(y(1)) - (sin((pi*sin(400*t*pi))/2)*sin(y(2)))/10 + 10*cos((pi*sin(400*t*pi))/2)*cos(y(1)) + (y(3)*cos(y(2))*sin(y(1)))/100 + (y(4)*cos(y(1))*sin(y(2)))/100 + (cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/10))/2000 - (y(4)*cos(y(1))*cos(y(2))*((sin((pi*sin(400*t*pi))/2)*sin(y(2)))/10 + y(3)*cos(y(1)) + 10*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - (y(3)*cos(y(2))*sin(y(1)))/100 - (y(4)*cos(y(1))*sin(y(2)))/100 - (cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/10))/2000 - (y(3)*sin(y(1))*sin(y(2))*(y(3)*cos(y(1)) - (sin((pi*sin(400*t*pi))/2)*sin(y(2)))/10 + 10*cos((pi*sin(400*t*pi))/2)*cos(y(1)) + (y(3)*cos(y(2))*sin(y(1)))/100 + (y(4)*cos(y(1))*sin(y(2)))/100 + (cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/10))/2000 + (y(3)*sin(y(1))*sin(y(2))*((sin((pi*sin(400*t*pi))/2)*sin(y(2)))/10 + y(3)*cos(y(1)) + 10*cos((pi*sin(400*t*pi))/2)*cos(y(1)) - (y(3)*cos(y(2))*sin(y(1)))/100 - (y(4)*cos(y(1))*sin(y(2)))/100 - (cos((pi*sin(400*t*pi))/2)*cos(y(2))*sin(y(1)))/10))/2000)/(((cos((pi*sin(400*t*pi))/2)*cos(y(2)))/100 - (sin((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/100)^2/10 + ((sin((pi*sin(400*t*pi))/2)*cos(y(2)))/100 + (cos((pi*sin(400*t*pi))/2)*sin(y(1))*sin(y(2)))/100)^2/10 + (cos(y(1))^2*sin(y(2))^2)/100000);
 t