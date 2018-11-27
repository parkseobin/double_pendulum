function runge_kutta(yin, l1, l2, h)
{
	var yout = [0, 0, 0, 0];
	var yt = [0, 0, 0, 0];
	var k = [[0, 0, 0, 0], 
	     [0, 0, 0, 0], 
		 [0, 0, 0, 0], 
		 [0, 0, 0, 0]];

	var i;

	var dydx = derivatives(yin, l1, l2);
	for(i=0; i<4; i++)
	{
		k[0][i] = h*dydx[i];
		yt[i] = yin[i] + 0.5*k[0][i];
	}

	dydx = derivatives(yt, l1, l2);
	for(i=0; i<4; i++)
	{
		k[1][i] = h*dydx[i];
		yt[i] = yin[i] + 0.5*k[1][i];
	}
	
	dydx = derivatives(yt, l1, l2);
	for(i=0; i<4; i++)
	{
		k[2][i] = h*dydx[i];
		yt[i] = yin[i] + k[2][i];
	}

	dydx = derivatives(yt, l1, l2);
	for(i=0; i<4; i++)
	{
		k[3][i] = h*dydx[i];
		yout[i] = yin[i] + k[0][i]/6. + k[1][i]/3. + k[2][i]/3. + k[3][i]/6.;
	}

	yout[0] = yout[0] % (Math.PI*2)
	yout[2] = yout[2] % (Math.PI*2)
	return yout;
};

function derivatives(yin, l1, l2){
	var m1 = 1.;
	var m2 = 1.;
	//l1 = 1.;
	//l2 = 1.;
	var G = 9.8;

	var dydx = [0, 0, 0, 0];
	dydx[0] = yin[1];

	var del = yin[2] - yin[0];
	var den1 = (m1+m2)*l1 - m2*l1*Math.cos(del)*Math.cos(del);
	dydx[1] = (m2*l1*yin[1]*yin[1]*Math.sin(del)*Math.cos(del)
			+ m2*G*Math.sin(yin[2])*Math.cos(del) 
			+ m2*l2*yin[3]*yin[3]*Math.sin(del)
			- (m1+m2)*G*Math.sin(yin[0]))/den1;

	dydx[2] = yin[3];

	var den2 = (l2/l1)*den1;
	dydx[3] = (-m2*l2*yin[3]*yin[3]*Math.sin(del)*Math.cos(del)
			+ (m1+m2)*G*Math.sin(yin[0])*Math.cos(del)
			- (m1+m2)*l1*yin[1]*yin[1]*Math.sin(del)
			- (m1+m2)*G*Math.sin(yin[2]))/den2;

	return dydx;
};

