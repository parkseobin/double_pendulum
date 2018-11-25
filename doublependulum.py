import numpy as np
import glfw
from OpenGL.GL import *
from OpenGL.GLU import *
from random import random

def render(th1, th2, camAng):

	glPolygonMode( GL_FRONT_AND_BACK, GL_LINE )
# enable depth test
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
	glEnable(GL_DEPTH_TEST)
	glLoadIdentity()

# orthogonal projection
	glOrtho(-1, 1, -1, 1, -3, 3)
# rotate camera position
	gluLookAt(.2*np.sin(camAng), .1,.2*np.cos(camAng), 0, 0, 0, 0, 1, 0)

# draw cooridnate
	glBegin(GL_LINES)
	glColor3ub(255, 0, 0)
	glVertex3fv(np.array([0.,0.,0.]))
	glVertex3fv(np.array([1.,0.,0.]))
	glColor3ub(0, 255, 0)
	glVertex3fv(np.array([0.,0.,0.]))
	glVertex3fv(np.array([0.,1.,0.]))
	glColor3ub(0, 0, 255)
	# aligning with pendulum path
	glVertex3fv(np.array([0.,0.,0.]))
	glVertex3fv(np.array([-1.,0.,0.]))
	glEnd()




	# start drawing pendulum
	
	#-------
	# LAYER1
	glPushMatrix()

	#-------
	# LAYER2
	glPushMatrix()

	glColor3ub(255, 0, 0)
	glScalef(.05, .05, .05)
	glColor3ub(255, 0, 0)

	# first ball
	drawSphere()
	glPopMatrix()
	glTranslatef(0.3*np.sin(th1), -0.3*np.cos(th1), 0)

	#-------
	# LAYER2
	glPushMatrix()
	glTranslatef(0.3*np.sin(th1), -0.3*np.cos(th1), 0)

	#-------
	# LAYER3
	glPushMatrix()
	glScalef(.05, .05, .05)
	glColor3ub(255, 0, 0)
	# second ball
	drawSphere()
	glPopMatrix()

	
	glTranslatef(0.3*np.sin(th1+th2), -0.3*np.cos(th1+th2), 0)

	#-------
	# LAYER3
	glPushMatrix()
	glTranslatef(0.3*np.sin(th1+th2), -0.3*np.cos(th1+th2), 0)
	glScalef(.05, .05, .05)
	glColor3ub(255, 0, 0)
	# third ball
	drawSphere()
	glPopMatrix()

	glRotatef(((th2+th1)*180/np.pi) % 360, 0, 0, 1)
	glColor3ub(0, 255, 0)
	# green cube
	glScalef(0.2, 2, 0.2)
	drawCube()

	glPopMatrix()


	glRotatef((th1*180/np.pi) % 360, 0, 0, 1)
	glColor3ub(0, 0, 255)
	# blue cube
	glScalef(0.2, 2, 0.2)
	drawCube()

	glPopMatrix()




# numerically compuing double pendulum with Runge-Kutta method
def runge_kutta(yin, h):
	yout = [0]*4
	yt = [0]*4
	k = np.zeros((4, 4))


	dydx = derivatives(yin)
	for i in range(4):
		k[0][i] = h*dydx[i]
		yt[i] = yin[i] + 0.5*k[0][i]


	dydx = derivatives(yt)
	for i in range(4):
		k[1][i] = h*dydx[i]
		yt[i] = yin[i] + 0.5*k[1][i]


	dydx = derivatives(yt)
	for i in range(4):
		k[2][i] = h*dydx[i]
		yt[i] = yin[i] + k[2][i]

	dydx = derivatives(yt)
	for i in range(4):
		k[3][i] = h*dydx[i]
		yout[i] = yin[i] + k[0][i]/6. + k[1][i]/3. + k[2][i]/3. + k[3][i]/6.

	return yout


def derivatives(yin):
	# masses and lengths of the pendulum, gravity constant
	M1 = 1.2
	M2 = 1.0
	L1 = 1.0
	L2 = 1.0
	G = 9.8
	
	dydx = [0]*4
	dydx[0] = yin[1]

	del_ = yin[2]-yin[0]
	den1 = (M1+M2)*L1 - M2*L1*np.cos(del_)*np.cos(del_);
	dydx[1] = (M2*L1*yin[1]*yin[1]*np.sin(del_)*np.cos(del_)
  		+ M2*G*np.sin(yin[2])*np.cos(del_) + M2*L2*yin[3]*yin[3]*np.sin(del_)
  		- (M1+M2)*G*np.sin(yin[0]))/den1	
	
	dydx[2] = yin[3]

	den2 = (L2/L1)*den1
	dydx[3] = (-M2*L2*yin[3]*yin[3]*np.sin(del_)*np.cos(del_)
  		+ (M1+M2)*G*np.sin(yin[0])*np.cos(del_) 
  		- (M1+M2)*L1*yin[1]*yin[1]*np.sin(del_)
  		- (M1+M2)*G*np.sin(yin[2]))/den2

	return dydx


def drawCube():
	glBegin(GL_QUADS)
	glVertex3f( 0.1, 0.1,-0.1)
	glVertex3f(-0.1, 0.1,-0.1)
	glVertex3f(-0.1, 0.1, 0.1)
	glVertex3f( 0.1, 0.1, 0.1)
	glVertex3f( 0.1,-0.1, 0.1)
	glVertex3f(-0.1,-0.1, 0.1)
	glVertex3f(-0.1,-0.1,-0.1)
	glVertex3f( 0.1,-0.1,-0.1)
	glVertex3f( 0.1, 0.1, 0.1)
	glVertex3f(-0.1, 0.1, 0.1)
	glVertex3f(-0.1,-0.1, 0.1)
	glVertex3f( 0.1,-0.1, 0.1)
	glVertex3f( 0.1,-0.1,-0.1)
	glVertex3f(-0.1,-0.1,-0.1)
	glVertex3f(-0.1, 0.1,-0.1)
	glVertex3f( 0.1, 0.1,-0.1)
	glVertex3f(-0.1, 0.1, 0.1)
	glVertex3f(-0.1, 0.1,-0.1)
	glVertex3f(-0.1,-0.1,-0.1)
	glVertex3f(-0.1,-0.1, 0.1)
	glVertex3f( 0.1, 0.1,-0.1)
	glVertex3f( 0.1, 0.1, 0.1)
	glVertex3f( 0.1,-0.1, 0.1)
	glVertex3f( 0.1,-0.1,-0.1)
	glEnd()


# draw a sphere of radius 1, centered at the origin.
# numLats: number of latitude segments (horizontal)
# numLongs: number of longitude segments (horizontal)
def drawSphere(numLats=12, numLongs=12):
	for i in range(0, numLats + 1):
		lat0 = np.pi * (-0.5 + float(float(i - 1) / float(numLats)))
		z0 = np.sin(lat0)
		zr0 = np.cos(lat0)
		lat1 = np.pi * (-0.5 + float(float(i) / float(numLats)))
		z1 = np.sin(lat1)
		zr1 = np.cos(lat1)
		# Use Quad strips to draw the sphere
		glBegin(GL_QUAD_STRIP)
		for j in range(0, numLongs + 1):
			lng = 2 * np.pi * float(float(j - 1) / float(numLongs))
			x = np.cos(lng)
			y = np.sin(lng)
			glVertex3f(x * zr0, y * zr0, z0)
			glVertex3f(x * zr1, y * zr1, z1)
		glEnd()


def key_callback(window, key, scancode, action, mods):
	global camAng
	if action==glfw.PRESS or action==glfw.REPEAT:
		if key==glfw.KEY_3:
			camAng += np.radians(5)
		if key==glfw.KEY_1:
			camAng -= np.radians(5)



def main():
	global camAng

	if not glfw.init():
		return
	window = glfw.create_window(700,700,"2014001303", None,None)
	glfw.set_key_callback(window, key_callback)
	if not window:
		glfw.terminate()
		return
	glfw.make_context_current(window)

	camAng = 0

	# randomly generates location of the pendulum
	th1 = np.pi * random()
	th2 = np.pi * random()
	w1 = 0
	w2 = 0
	while not glfw.window_should_close(window):
		glfw.poll_events()

		render(th1, th2, camAng)
		(th1, w1, th2, w2) = runge_kutta([th1, w1, th2, w2], 0.003)
		glfw.swap_buffers(window)
		
	glfw.terminate()


if __name__ == "__main__":
	main()
