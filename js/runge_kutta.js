function runge_kutta(yin, h)
{
	yout = [0, 0, 0, 0];
	yt = [0, 0, 0, 0];
	k = [[0, 0, 0, 0], 
	     [0, 0, 0, 0], 
		 [0, 0, 0, 0], 
		 [0, 0, 0, 0]];

	var i;

	dydx = derivatives(yin);
	for(i=0; i<4; i++)
	{
		k[0][i] = h*dydx[i];
		yt[i] = yin[i] + 0.5*k[0][i];
	}

	dydx = derivatives(yt);
	for(i=0; i<4; i++)
	{
		k[1][i] = h*dydx[i];
		yt[i] = yin[i] + 0.5*k[1][i];
	}
	
	dydx = derivatives(yt);
	for(i=0; i<4; i++)
	{
		k[2][i] = h*dydx[i];
		yt[i] = yin[i] + k[2][i];
	}

	dydx = derivatives(yt);
	for(i=0; i<4; i++)
	{
		k[3][i] = h*dydx[i];
		yout[i] = yin[i] + k[0][i]/6. + k[1][i]/3. + k[2][i]/3. + k[3][i]/6.;
	}
	return yout;
};

function derivatives(yin){
	m1 = 1.;
	m2 = 1.;
	l1 = 1.;
	l2 = 1.;
	G = 9.8;

	dydx = [0, 0, 0, 0];
	dydx[0] = yin[1];

	del = yin[2] - yin[0];
	den1 = (m1+m2)*l1 - m2*l1*Math.cos(del)*Math.cos(del);
	dydx[1] = (m2*l1*yin[1]*yin[1]*Math.sin(del)*Math.cos(del)
			+ m2*G*Math.sin(yin[2])*Math.cos(del) 
			+ m2*l2*yin[3]*yin[3]*Math.sin(del)
			- (m1+m2)*G*Math.sin(yin[0]))/den1;

	dydx[2] = yin[3];

	den2 = (l2/l1)*den1;
	dydx[3] = (-m2*l2*yin[3]*yin[3]*Math.sin(del)*Math.cos(del)
			+ (m1+m2)*G*Math.sin(yin[0])*Math.cos(del)
			- (m1+m2)*l1*yin[1]*yin[1]*Math.sin(del)
			- (m1+m2)*G*Math.sin(yin[2]))/den2;

	return dydx;
};

