function dy = myODE_twist(t,y)
 dy(1,1) =y(3);
 dy(2,1) =y(4);
 dy(3,1) =-(((cos((pi*sin(80*t*pi))/2)*sin(y(1)) - (cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/100)*(10*y(3)*cos(y(1)) + (y(3)*cos(y(2))*sin(y(1)))/10 + (y(4)*cos(y(1))*sin(y(2)))/10 + y(3)^2*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - (y(4)^2*sin((pi*sin(80*t*pi))/2)*sin(y(2)))/100 - 3200*pi^3*sin(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1)) - 32*pi^3*sin(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(2)) + 1600*pi^4*cos(80*t*pi)^2*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - 16*pi^4*cos(80*t*pi)^2*sin((pi*sin(80*t*pi))/2)*sin(y(2)) + (y(3)^2*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)^2*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 - 32*pi^3*sin(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)) + 16*pi^4*cos(80*t*pi)^2*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)) + (4*pi^2*y(4)*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(2)))/5 - 80*pi^2*y(3)*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(1)) + (y(3)*y(4)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/50 + (4*pi^2*y(3)*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/5 - (4*pi^2*y(4)*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/5))/4 - ((10*cos(y(1)) - 40*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(1)))*(10*sin(y(1)) + 40*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1))))/2 - ((cos((pi*sin(80*t*pi))/2)*sin(y(1)) + (cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/100)*((y(3)*cos(y(2))*sin(y(1)))/10 - 10*y(3)*cos(y(1)) + (y(4)*cos(y(1))*sin(y(2)))/10 - y(3)^2*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - (y(4)^2*sin((pi*sin(80*t*pi))/2)*sin(y(2)))/100 + 3200*pi^3*sin(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1)) - 32*pi^3*sin(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(2)) - 1600*pi^4*cos(80*t*pi)^2*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - 16*pi^4*cos(80*t*pi)^2*sin((pi*sin(80*t*pi))/2)*sin(y(2)) + (y(3)^2*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)^2*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 - 32*pi^3*sin(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)) + 16*pi^4*cos(80*t*pi)^2*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)) + (4*pi^2*y(4)*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(2)))/5 + 80*pi^2*y(3)*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(1)) + (y(3)*y(4)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/50 + (4*pi^2*y(3)*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/5 - (4*pi^2*y(4)*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/5))/4 + (5493695805276885*cos(y(1))*sin(y(1)))/1073741824 + 144*pi*y(3) - ((cos(y(1)) - (cos(y(2))*sin(y(1)))/100)*(y(3)^2*sin(y(1)) + 10*y(3)*cos((pi*sin(80*t*pi))/2)*sin(y(1)) - (y(4)*sin((pi*sin(80*t*pi))/2)*cos(y(2)))/10 + (y(3)^2*cos(y(1))*cos(y(2)))/100 + (y(4)^2*cos(y(1))*cos(y(2)))/100 - (y(3)*y(4)*sin(y(1))*sin(y(2)))/50 + 400*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1)) - 4*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(2)) + (y(3)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/10 - (y(4)*cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/10 - 4*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1))))/4 - ((cos(y(1)) + (cos(y(2))*sin(y(1)))/100)*(y(3)^2*sin(y(1)) + 10*y(3)*cos((pi*sin(80*t*pi))/2)*sin(y(1)) + (y(4)*sin((pi*sin(80*t*pi))/2)*cos(y(2)))/10 - (y(3)^2*cos(y(1))*cos(y(2)))/100 - (y(4)^2*cos(y(1))*cos(y(2)))/100 + (y(3)*y(4)*sin(y(1))*sin(y(2)))/50 + 400*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1)) + 4*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(2)) - (y(3)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/10 + (y(4)*cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/10 + 4*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1))))/4 + ((sin((pi*sin(80*t*pi))/2)*sin(y(1)) - (sin((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/100)*(y(3)^2*sin((pi*sin(80*t*pi))/2)*cos(y(1)) + (y(4)^2*cos((pi*sin(80*t*pi))/2)*sin(y(2)))/100 + 3200*pi^3*sin(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - 32*pi^3*sin(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(2)) + 1600*pi^4*cos(80*t*pi)^2*sin((pi*sin(80*t*pi))/2)*cos(y(1)) + 16*pi^4*cos(80*t*pi)^2*cos((pi*sin(80*t*pi))/2)*sin(y(2)) + (y(3)^2*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)^2*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + 32*pi^3*sin(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)) + 16*pi^4*cos(80*t*pi)^2*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)) + 80*pi^2*y(3)*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(1)) + (4*pi^2*y(4)*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2)))/5 + (y(3)*y(4)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/50 - (4*pi^2*y(3)*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/5 + (4*pi^2*y(4)*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/5))/4 - ((sin((pi*sin(80*t*pi))/2)*sin(y(1)) + (sin((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/100)*((y(4)^2*cos((pi*sin(80*t*pi))/2)*sin(y(2)))/100 - y(3)^2*sin((pi*sin(80*t*pi))/2)*cos(y(1)) - 3200*pi^3*sin(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - 32*pi^3*sin(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(2)) - 1600*pi^4*cos(80*t*pi)^2*sin((pi*sin(80*t*pi))/2)*cos(y(1)) + 16*pi^4*cos(80*t*pi)^2*cos((pi*sin(80*t*pi))/2)*sin(y(2)) + (y(3)^2*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)^2*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + 32*pi^3*sin(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)) + 16*pi^4*cos(80*t*pi)^2*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)) - 80*pi^2*y(3)*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(1)) + (4*pi^2*y(4)*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2)))/5 + (y(3)*y(4)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/50 - (4*pi^2*y(3)*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/5 + (4*pi^2*y(4)*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/5))/4 - (cos(y(1))*(y(3)^2*sin(y(1)) + 10*y(3)*cos((pi*sin(80*t*pi))/2)*sin(y(1)) + 400*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1))))/2 + ((y(3)*sin((pi*sin(80*t*pi))/2)*cos(y(1)) + 40*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(1)) - (y(3)*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 - (y(4)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/5)*((y(4)*cos((pi*sin(80*t*pi))/2)*cos(y(2)))/100 + y(3)*sin((pi*sin(80*t*pi))/2)*sin(y(1)) - 40*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(2)))/5 + (y(3)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/100 - (y(4)*sin((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/100 + (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/5))/4 - ((y(3)*sin((pi*sin(80*t*pi))/2)*cos(y(1)) + 40*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(1)) + (y(3)*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/100 - (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/5)*((y(4)*cos((pi*sin(80*t*pi))/2)*cos(y(2)))/100 - y(3)*sin((pi*sin(80*t*pi))/2)*sin(y(1)) + 40*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(2)))/5 + (y(3)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/100 - (y(4)*sin((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/100 + (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/5))/4 + ((y(3)*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - 40*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(1)) + (y(3)*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/5)*(10*sin(y(1)) - (cos(y(1))*cos(y(2)))/10 + y(3)*cos((pi*sin(80*t*pi))/2)*sin(y(1)) + (y(4)*sin((pi*sin(80*t*pi))/2)*cos(y(2)))/100 + 40*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1)) + (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(2)))/5 - (y(3)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/100 + (y(4)*cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/100 + (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/5))/4 - ((40*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(1)) - y(3)*cos((pi*sin(80*t*pi))/2)*cos(y(1)) + (y(3)*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/5)*(10*sin(y(1)) + (cos(y(1))*cos(y(2)))/10 + y(3)*cos((pi*sin(80*t*pi))/2)*sin(y(1)) - (y(4)*sin((pi*sin(80*t*pi))/2)*cos(y(2)))/100 + 40*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1)) - (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(2)))/5 + (y(3)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/100 - (y(4)*cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/100 - (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/5))/4 - (((cos(y(2))*sin(y(1)))/10 - 10*cos(y(1)) + 40*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(1)) + (y(4)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/5)*((y(4)*sin((pi*sin(80*t*pi))/2)*cos(y(2)))/100 - (cos(y(1))*cos(y(2)))/10 - 10*sin(y(1)) - 40*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1)) + (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(2)))/5 + (y(4)*cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/100 + (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/5))/4 - ((10*cos(y(1)) + (cos(y(2))*sin(y(1)))/10 - 40*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(1)) + (y(4)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/5)*(10*sin(y(1)) - (cos(y(1))*cos(y(2)))/10 + (y(4)*sin((pi*sin(80*t*pi))/2)*cos(y(2)))/100 + 40*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1)) + (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(2)))/5 + (y(4)*cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/100 + (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/5))/4 + ((40*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(1)) - (y(4)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/5)*(40*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - (y(4)*cos((pi*sin(80*t*pi))/2)*cos(y(2)))/100 + (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(2)))/5 + (y(4)*sin((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/100 - (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/5))/4 + ((40*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(1)) + (y(4)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/100 - (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/5)*((y(4)*cos((pi*sin(80*t*pi))/2)*cos(y(2)))/100 + 40*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(2)))/5 - (y(4)*sin((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/100 + (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/5))/4 - ((y(3)*sin(y(1)) - (y(3)*cos(y(1))*cos(y(2)))/100 + (y(4)*sin(y(1))*sin(y(2)))/100)*(y(3)*cos(y(1)) - (sin((pi*sin(80*t*pi))/2)*sin(y(2)))/10 + 10*cos((pi*sin(80*t*pi))/2)*cos(y(1)) + (y(3)*cos(y(2))*sin(y(1)))/100 + (y(4)*cos(y(1))*sin(y(2)))/100 + (cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/10))/4 - ((y(3)*sin(y(1)) + (y(3)*cos(y(1))*cos(y(2)))/100 - (y(4)*sin(y(1))*sin(y(2)))/100)*((sin((pi*sin(80*t*pi))/2)*sin(y(2)))/10 + y(3)*cos(y(1)) + 10*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - (y(3)*cos(y(2))*sin(y(1)))/100 - (y(4)*cos(y(1))*sin(y(2)))/100 - (cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/10))/4 + ((10*cos((pi*sin(80*t*pi))/2)*sin(y(1)) - (y(4)*sin(y(1))*sin(y(2)))/100 + (cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/10)*((sin((pi*sin(80*t*pi))/2)*sin(y(2)))/10 + 10*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - (y(4)*cos(y(1))*sin(y(2)))/100 - (cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/10))/4 + ((10*cos((pi*sin(80*t*pi))/2)*sin(y(1)) + (y(4)*sin(y(1))*sin(y(2)))/100 - (cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/10)*(10*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - (sin((pi*sin(80*t*pi))/2)*sin(y(2)))/10 + (y(4)*cos(y(1))*sin(y(2)))/100 + (cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/10))/4 - (5493642118185685*cos((pi*sin(80*t*pi))/2)^2*cos(y(1))*sin(y(1)))/1073741824 + (cos((pi*sin(80*t*pi))/2)*sin(y(1))*(10*y(3)*cos(y(1)) + y(3)^2*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - 3200*pi^3*sin(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1)) + 1600*pi^4*cos(80*t*pi)^2*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - 80*pi^2*y(3)*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(1))))/2 + (sin((pi*sin(80*t*pi))/2)*sin(y(1))*(y(3)^2*sin((pi*sin(80*t*pi))/2)*cos(y(1)) + 3200*pi^3*sin(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1)) + 1600*pi^4*cos(80*t*pi)^2*sin((pi*sin(80*t*pi))/2)*cos(y(1)) + 80*pi^2*y(3)*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(1))))/2 - (y(3)*sin(y(1))*(y(3)*cos(y(1)) + 10*cos((pi*sin(80*t*pi))/2)*cos(y(1))))/2 + (y(3)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*(y(3)*sin((pi*sin(80*t*pi))/2)*sin(y(1)) - 40*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1))))/2 + (y(3)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*(10*sin(y(1)) + y(3)*cos((pi*sin(80*t*pi))/2)*sin(y(1)) + 40*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1))))/2 + 20*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(1))*(y(3)*sin((pi*sin(80*t*pi))/2)*sin(y(1)) - 40*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1))) - 20*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(1))*(10*sin(y(1)) + y(3)*cos((pi*sin(80*t*pi))/2)*sin(y(1)) + 40*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1))) + 800*pi^4*cos(80*t*pi)^2*cos((pi*sin(80*t*pi))/2)^2*cos(y(1))*sin(y(1)))/(cos(y(1))^2/2 + (cos((pi*sin(80*t*pi))/2)*sin(y(1)) - (cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/100)^2/4 + (cos((pi*sin(80*t*pi))/2)*sin(y(1)) + (cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/100)^2/4 + (sin((pi*sin(80*t*pi))/2)*sin(y(1)) - (sin((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/100)^2/4 + (sin((pi*sin(80*t*pi))/2)*sin(y(1)) + (sin((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/100)^2/4 + (cos((pi*sin(80*t*pi))/2)^2*sin(y(1))^2)/2 + (cos(y(1)) - (cos(y(2))*sin(y(1)))/100)^2/4 + (cos(y(1)) + (cos(y(2))*sin(y(1)))/100)^2/4 + (sin((pi*sin(80*t*pi))/2)^2*sin(y(1))^2)/2);
 dy(4,1) =-((5493695805276885*y(2))/2147483648 + 144*pi*y(4) + (((sin((pi*sin(80*t*pi))/2)*cos(y(2)))/100 + (cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/100)*(10*y(3)*cos(y(1)) + (y(3)*cos(y(2))*sin(y(1)))/10 + (y(4)*cos(y(1))*sin(y(2)))/10 + y(3)^2*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - (y(4)^2*sin((pi*sin(80*t*pi))/2)*sin(y(2)))/100 - 3200*pi^3*sin(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1)) - 32*pi^3*sin(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(2)) + 1600*pi^4*cos(80*t*pi)^2*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - 16*pi^4*cos(80*t*pi)^2*sin((pi*sin(80*t*pi))/2)*sin(y(2)) + (y(3)^2*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)^2*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 - 32*pi^3*sin(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)) + 16*pi^4*cos(80*t*pi)^2*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)) + (4*pi^2*y(4)*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(2)))/5 - 80*pi^2*y(3)*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(1)) + (y(3)*y(4)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/50 + (4*pi^2*y(3)*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/5 - (4*pi^2*y(4)*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/5))/4 + (((sin((pi*sin(80*t*pi))/2)*cos(y(2)))/100 + (cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/100)*((y(3)*cos(y(2))*sin(y(1)))/10 - 10*y(3)*cos(y(1)) + (y(4)*cos(y(1))*sin(y(2)))/10 - y(3)^2*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - (y(4)^2*sin((pi*sin(80*t*pi))/2)*sin(y(2)))/100 + 3200*pi^3*sin(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1)) - 32*pi^3*sin(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(2)) - 1600*pi^4*cos(80*t*pi)^2*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - 16*pi^4*cos(80*t*pi)^2*sin((pi*sin(80*t*pi))/2)*sin(y(2)) + (y(3)^2*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)^2*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 - 32*pi^3*sin(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)) + 16*pi^4*cos(80*t*pi)^2*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)) + (4*pi^2*y(4)*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(2)))/5 + 80*pi^2*y(3)*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(1)) + (y(3)*y(4)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/50 + (4*pi^2*y(3)*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/5 - (4*pi^2*y(4)*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/5))/4 - (((cos((pi*sin(80*t*pi))/2)*cos(y(2)))/100 - (sin((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/100)*(y(3)^2*sin((pi*sin(80*t*pi))/2)*cos(y(1)) + (y(4)^2*cos((pi*sin(80*t*pi))/2)*sin(y(2)))/100 + 3200*pi^3*sin(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - 32*pi^3*sin(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(2)) + 1600*pi^4*cos(80*t*pi)^2*sin((pi*sin(80*t*pi))/2)*cos(y(1)) + 16*pi^4*cos(80*t*pi)^2*cos((pi*sin(80*t*pi))/2)*sin(y(2)) + (y(3)^2*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)^2*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + 32*pi^3*sin(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)) + 16*pi^4*cos(80*t*pi)^2*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)) + 80*pi^2*y(3)*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(1)) + (4*pi^2*y(4)*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2)))/5 + (y(3)*y(4)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/50 - (4*pi^2*y(3)*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/5 + (4*pi^2*y(4)*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/5))/4 - (((cos((pi*sin(80*t*pi))/2)*cos(y(2)))/100 - (sin((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/100)*((y(4)^2*cos((pi*sin(80*t*pi))/2)*sin(y(2)))/100 - y(3)^2*sin((pi*sin(80*t*pi))/2)*cos(y(1)) - 3200*pi^3*sin(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - 32*pi^3*sin(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(2)) - 1600*pi^4*cos(80*t*pi)^2*sin((pi*sin(80*t*pi))/2)*cos(y(1)) + 16*pi^4*cos(80*t*pi)^2*cos((pi*sin(80*t*pi))/2)*sin(y(2)) + (y(3)^2*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (y(4)^2*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + 32*pi^3*sin(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)) + 16*pi^4*cos(80*t*pi)^2*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)) - 80*pi^2*y(3)*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(1)) + (4*pi^2*y(4)*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2)))/5 + (y(3)*y(4)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/50 - (4*pi^2*y(3)*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/5 + (4*pi^2*y(4)*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/5))/4 - (((cos(y(1))*sin(y(2)))/10 + (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(2)))/5 + (y(3)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/100 - (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/5)*(10*sin(y(1)) - (cos(y(1))*cos(y(2)))/10 + y(3)*cos((pi*sin(80*t*pi))/2)*sin(y(1)) + 40*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1)) + (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(2)))/5 - (y(3)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/100 + (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/5))/4 + (((cos(y(1))*sin(y(2)))/10 + (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(2)))/5 + (y(3)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/100 - (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/5)*(10*sin(y(1)) + (cos(y(1))*cos(y(2)))/10 + y(3)*cos((pi*sin(80*t*pi))/2)*sin(y(1)) + 40*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1)) - (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(2)))/5 + (y(3)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/100 - (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/5))/4 - (((y(4)*cos((pi*sin(80*t*pi))/2)*sin(y(2)))/100 + (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2)))/5 + (y(3)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + (y(4)*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/5)*((y(4)*cos((pi*sin(80*t*pi))/2)*cos(y(2)))/100 + y(3)*sin((pi*sin(80*t*pi))/2)*sin(y(1)) - 40*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(2)))/5 + (y(3)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/100 - (y(4)*sin((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/100 + (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/5))/4 - (((y(4)*cos((pi*sin(80*t*pi))/2)*sin(y(2)))/100 + (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2)))/5 + (y(3)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + (y(4)*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 + (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/5)*((y(4)*cos((pi*sin(80*t*pi))/2)*cos(y(2)))/100 - y(3)*sin((pi*sin(80*t*pi))/2)*sin(y(1)) + 40*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(2)))/5 + (y(3)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/100 - (y(4)*sin((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/100 + (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/5))/4 + (((sin((pi*sin(80*t*pi))/2)*cos(y(2)))/10 + (y(3)*sin(y(1))*sin(y(2)))/100 + (cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/10)*(y(3)*cos(y(1)) - (sin((pi*sin(80*t*pi))/2)*sin(y(2)))/10 + 10*cos((pi*sin(80*t*pi))/2)*cos(y(1)) + (y(3)*cos(y(2))*sin(y(1)))/100 + (cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/10))/4 - (((sin((pi*sin(80*t*pi))/2)*cos(y(2)))/10 + (y(3)*sin(y(1))*sin(y(2)))/100 + (cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/10)*((sin((pi*sin(80*t*pi))/2)*sin(y(2)))/10 + y(3)*cos(y(1)) + 10*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - (y(3)*cos(y(2))*sin(y(1)))/100 - (cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/10))/4 + (((2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(2)))/5 - (y(4)*sin((pi*sin(80*t*pi))/2)*sin(y(2)))/100 + (y(3)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + (y(4)*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 - (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/5)*(10*sin(y(1)) - (cos(y(1))*cos(y(2)))/10 + y(3)*cos((pi*sin(80*t*pi))/2)*sin(y(1)) + (y(4)*sin((pi*sin(80*t*pi))/2)*cos(y(2)))/100 + 40*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1)) + (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(2)))/5 - (y(3)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/100 + (y(4)*cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/100 + (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/5))/4 - (((2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(2)))/5 - (y(4)*sin((pi*sin(80*t*pi))/2)*sin(y(2)))/100 + (y(3)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + (y(4)*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/100 - (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/5)*(10*sin(y(1)) + (cos(y(1))*cos(y(2)))/10 + y(3)*cos((pi*sin(80*t*pi))/2)*sin(y(1)) - (y(4)*sin((pi*sin(80*t*pi))/2)*cos(y(2)))/100 + 40*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1)) - (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(2)))/5 + (y(3)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/100 - (y(4)*cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/100 - (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/5))/4 + (((2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2)))/5 + (y(3)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/5)*(y(3)*sin((pi*sin(80*t*pi))/2)*sin(y(1)) - 40*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(2)))/5 + (y(3)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/100 + (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/5))/4 + (((2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2)))/5 + (y(3)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*sin(y(2)))/100 + (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/5)*(40*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - y(3)*sin((pi*sin(80*t*pi))/2)*sin(y(1)) - (2*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*sin(y(2)))/5 + (y(3)*sin((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/100 + (2*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/5))/4 + (cos(y(1))*sin(y(2))*(y(3)^2*sin(y(1)) + 10*y(3)*cos((pi*sin(80*t*pi))/2)*sin(y(1)) - (y(4)*sin((pi*sin(80*t*pi))/2)*cos(y(2)))/10 + (y(3)^2*cos(y(1))*cos(y(2)))/100 + (y(4)^2*cos(y(1))*cos(y(2)))/100 - (y(3)*y(4)*sin(y(1))*sin(y(2)))/50 + 400*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1)) - 4*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(2)) + (y(3)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/10 - (y(4)*cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/10 - 4*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1))))/400 - (cos(y(1))*sin(y(2))*(y(3)^2*sin(y(1)) + 10*y(3)*cos((pi*sin(80*t*pi))/2)*sin(y(1)) + (y(4)*sin((pi*sin(80*t*pi))/2)*cos(y(2)))/10 - (y(3)^2*cos(y(1))*cos(y(2)))/100 - (y(4)^2*cos(y(1))*cos(y(2)))/100 + (y(3)*y(4)*sin(y(1))*sin(y(2)))/50 + 400*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(1)) + 4*pi^2*cos(80*t*pi)*cos((pi*sin(80*t*pi))/2)*sin(y(2)) - (y(3)*cos((pi*sin(80*t*pi))/2)*cos(y(1))*cos(y(2)))/10 + (y(4)*cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/10 + 4*pi^2*cos(80*t*pi)*sin((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1))))/400 + (y(4)*cos(y(1))*cos(y(2))*(y(3)*cos(y(1)) - (sin((pi*sin(80*t*pi))/2)*sin(y(2)))/10 + 10*cos((pi*sin(80*t*pi))/2)*cos(y(1)) + (y(3)*cos(y(2))*sin(y(1)))/100 + (y(4)*cos(y(1))*sin(y(2)))/100 + (cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/10))/400 - (y(4)*cos(y(1))*cos(y(2))*((sin((pi*sin(80*t*pi))/2)*sin(y(2)))/10 + y(3)*cos(y(1)) + 10*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - (y(3)*cos(y(2))*sin(y(1)))/100 - (y(4)*cos(y(1))*sin(y(2)))/100 - (cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/10))/400 - (y(3)*sin(y(1))*sin(y(2))*(y(3)*cos(y(1)) - (sin((pi*sin(80*t*pi))/2)*sin(y(2)))/10 + 10*cos((pi*sin(80*t*pi))/2)*cos(y(1)) + (y(3)*cos(y(2))*sin(y(1)))/100 + (y(4)*cos(y(1))*sin(y(2)))/100 + (cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/10))/400 + (y(3)*sin(y(1))*sin(y(2))*((sin((pi*sin(80*t*pi))/2)*sin(y(2)))/10 + y(3)*cos(y(1)) + 10*cos((pi*sin(80*t*pi))/2)*cos(y(1)) - (y(3)*cos(y(2))*sin(y(1)))/100 - (y(4)*cos(y(1))*sin(y(2)))/100 - (cos((pi*sin(80*t*pi))/2)*cos(y(2))*sin(y(1)))/10))/400)/(((cos((pi*sin(80*t*pi))/2)*cos(y(2)))/100 - (sin((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/100)^2/2 + ((sin((pi*sin(80*t*pi))/2)*cos(y(2)))/100 + (cos((pi*sin(80*t*pi))/2)*sin(y(1))*sin(y(2)))/100)^2/2 + (cos(y(1))^2*sin(y(2))^2)/20000);
 t