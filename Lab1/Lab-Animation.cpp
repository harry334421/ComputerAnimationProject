#include <assert.h>
#include <math.h>

// glut
#include <GL/glut.h>

//================================
// global variables
//================================
// screen size
int g_screenWidth = 0;
int g_screenHeight = 0;

//number of Points for Spline
static int points = 0; //points index from 0
static int number = 8; //number of points

//time variables Q(t)
static GLfloat t = 0;

//M matrix
static GLfloat M[16] = { 0 };

//Catmul-Rom Spline M Marix
static GLfloat CRSplineM[16] = { -0.5f, 1.0f, -0.5f, 0.0f,	  // Column 1
								1.5f, -2.5f, 0.0f, 1.0f,      // Column 2
								-1.5f, 2.0f, 0.5f, 0.0f,      // Column 3
								0.5f, -0.5f, 0.0f, 0.0f };    // Column 4

//B Spline M Marix
static GLfloat BSplineM[16] = { -1.0f / 6.0f, 3.0f / 6.0f, -3.0f / 6.0f, 1.0f / 6.0f,  // Column 1
								3.0f / 6.0f, -6.0f / 6.0f, 0.0f / 6.0f, 4.0f / 6.0f,  // Column 2
								-3.0f / 6.0f, 3.0f / 6.0f, 3.0f / 6.0f, 1.0f / 6.0f,  // Column 3
								1.0f / 6.0f, 0.0f / 6.0f, 0.0f / 6.0f, 0.0f / 6.0f };// Column 4

//w_q, x_q, y_q, z_q, x, y, z
static GLfloat points_quater[8][7] = { { 1, 0, 0, 0,-5, 0, -5 },           //point 1
											{ 0, 1, 0, 0, -2, 2, -10 },    //point 2
											{ 0, 0, 1, 0, -1, 1, -15 },    //point 3
											{ 0, 0, 0, 1, 0, -5, -20 },    //point 4
											{ 0, 0, 1, 0, 1, 1, -15 },     //point 5
											{ 0, 1, 0, 0, 2, 2, -10 },     //point 6
											{ 1, 0, 0, 0, 5, 0, -5 },      //point 7
											{ 0, 1, 0, 0, 1, -3,-10} };   //point 8

//x_angle, y_angle, z_angle, x, y, z
static GLfloat points_euler[8][6] = { { 90, 0, 45, -5, 0, -5 },             //point 1
											{ 70, 20, 65, -2, 2, -10 },	    //point 2
											{ 50, 40, 85, -1, 1, -15 },	    //point 3
											{ 30, 60, 105, 0, -5, -20 },	//point 4
											{ 50, 40, 85, 1, 1, -15 },   	//point 5
											{ 70, 20, 65, 3, 3, -10 },		//point 6
											{ 90, 0, 45, 5, 0, -5 },    	//point 7
											{ 100, 20, 65, 1, -3, -10 } };  //point 8

//Q(t)=T*M*G
GLfloat mix(GLfloat T[4], GLfloat M[16], GLfloat G[4]) {
	GLfloat R[4] = { 0 };
	R[0] = T[0] * M[0] + T[1] * M[1] + T[2] * M[2] + T[3] * M[3];
	R[1] = T[0] * M[4] + T[1] * M[5] + T[2] * M[6] + T[3] * M[7];
	R[2] = T[0] * M[8] + T[1] * M[9] + T[2] * M[10] + T[3] * M[11];
	R[3] = T[0] * M[12] + T[1] * M[13] + T[2] * M[14] + T[3] * M[15];

	GLfloat Q = R[0] * G[0] + R[1] * G[1] + R[2] * G[2] + R[3] * G[3];

	return Q;
}

void Normalization(GLfloat N_tempM[8]) {
	GLfloat squar_quaterion = N_tempM[0] * N_tempM[0] + N_tempM[1] * N_tempM[1] + N_tempM[2] * N_tempM[2] + N_tempM[3] * N_tempM[3];
	if (squar_quaterion != 0) // ensure it's not divided by 0
	{
		GLfloat base_quaternion = sqrt(squar_quaterion);
		N_tempM[0] = N_tempM[0] / base_quaternion;
		N_tempM[1] = N_tempM[1] / base_quaternion;
		N_tempM[2] = N_tempM[2] / base_quaternion;
		N_tempM[3] = N_tempM[3] / base_quaternion;
	}
}

void QuaternionRoatationM(GLfloat Q_M[7], GLfloat R[16])
{
	GLfloat w = Q_M[0];
	GLfloat x = Q_M[1];
	GLfloat y = Q_M[2];
	GLfloat z = Q_M[3];
	R[0] = 1.0f - 2.0f * y * y - 2.0f * z * z;    //column1 row1
	R[1] = 2.0f * x * y + 2.0f * w * z;           //....... row2
	R[2] = 2.0f * x * z - 2.0f * w * y;		      //....... row3
	R[3] = 0.0f;					              //....... row4
	R[4] = 2.0f * x * y - 2.0f * w * z;		      //column2 row1
	R[5] = 1.0f - 2.0f * x * x - 2.0f * z * z;    //....... row2
	R[6] = 2.0f * y * z + 2.0f * w * x;		      //....... row3
	R[7] = 0.0f;					              //....... row4
	R[8] = 2.0f * x * z + 2.0f * w * y;		      //column3 row1
	R[9] = 2.0f * y * z - 2.0f * w * x;		      //....... row2
	R[10] = 1.0f - 2.0f * x * x - 2.0f * y * y;   //....... row3
	R[11] = 0.0f;                                 //....... row4
	R[12] = Q_M[4];				                  //column4 row1
	R[13] = Q_M[5];			                      //....... row2
	R[14] = Q_M[6];			                      //....... row3
	R[15] = 1.0f;					              //....... row4
}

void Euler2Quaternion(GLfloat E_M[8])
{
	GLfloat a = E_M[0] / 2;
	GLfloat b = E_M[1] / 2;
	GLfloat c = E_M[2] / 2;

	// // put each value 1 position next to original, since 3 Euler Angle turned to 4 Quaternion
	E_M[6] = E_M[5];
	E_M[5] = E_M[4];
	E_M[4] = E_M[3];
	E_M[0] = cos(c) * cos(b) * cos(c) + sin(c) * sin(b) * sin(a); //w
	E_M[1] = sin(c) * cos(b) * cos(c) - cos(c) * sin(b) * sin(a); //x
	E_M[2] = cos(c) * sin(b) * cos(c) + sin(c) * cos(b) * sin(a); //y
	E_M[3] = cos(c) * cos(b) * sin(c) - sin(c) * sin(b) * cos(a); //z
}

//Quaternion Interpolating
void q_interpolate(GLfloat p_quaternion[6][7], GLfloat SplineM[16]) {
	//T matrix
	GLfloat TMatrix_q[4] = { t * t * t, t * t, t, 1 };

	GLfloat tempM[7];

	for (int i = 0; i < 7; i++)
	{
		// changed by timer
		GLfloat GMatrix_q[4] = { p_quaternion[points][i],
								 p_quaternion[(points + 1)][i],
								 p_quaternion[(points + 2)][i],
								 p_quaternion[(points + 3)][i] };

		tempM[i] = mix(TMatrix_q, SplineM, GMatrix_q);
	}
	Normalization(tempM);
	QuaternionRoatationM(tempM, M);
}

//Euler Angle Interpolating
void e_interpolate(GLfloat p_euler[7][6], GLfloat SplineM[16])
{
	//T matrix T
	GLfloat TMatrix_e[4] = { t * t * t, t * t, t, 1 };

	GLfloat tempM[7];

	for (int i = 0; i < 7; i++)
	{
		// changed by timer
		GLfloat GMatrix_e[4] = { p_euler[points][i],
								 p_euler[(points + 1)][i],
								 p_euler[(points + 2)][i],
								 p_euler[(points + 3)][i] };

		tempM[i] = mix(TMatrix_e, SplineM, GMatrix_e);
	}
	Euler2Quaternion(tempM);
	Normalization(tempM);
	QuaternionRoatationM(tempM, M);
}

void teapot()
{
	/*change for different animation effect*/

	//q_interpolate(points_quater, CRSplineM);
	//q_interpolate(points_quater, BSplineM);
	//e_interpolate(points_euler, CRSplineM);
	e_interpolate(points_euler, BSplineM);
	glLoadMatrixf(M);

	// render objects
	glutSolidTeapot(1.0);

}

//================================
// render
//================================
void render(void) {
	// clear buffer
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClearDepth(1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// render state
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);

	// enable lighting
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	// light source attributes
	GLfloat LightAmbient[] = { 0.4f, 0.4f, 0.4f, 1.0f };
	GLfloat LightDiffuse[] = { 0.3f, 0.3f, 0.3f, 1.0f };
	GLfloat LightSpecular[] = { 0.4f, 0.4f, 0.4f, 1.0f };
	GLfloat LightPosition[] = { 5.0f, 5.0f, 5.0f, 1.0f };

	glLightfv(GL_LIGHT0, GL_AMBIENT, LightAmbient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, LightDiffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, LightSpecular);
	glLightfv(GL_LIGHT0, GL_POSITION, LightPosition);

	// surface material attributes
	GLfloat material_Ka[] = { 0.11f, 0.06f, 0.11f, 1.0f };
	GLfloat material_Kd[] = { 0.43f, 0.47f, 0.54f, 1.0f };
	GLfloat material_Ks[] = { 0.33f, 0.33f, 0.52f, 1.0f };
	GLfloat material_Ke[] = { 0.1f , 0.0f , 0.1f , 1.0f };
	GLfloat material_Se = 10;

	glMaterialfv(GL_FRONT, GL_AMBIENT, material_Ka);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, material_Kd);
	glMaterialfv(GL_FRONT, GL_SPECULAR, material_Ks);
	glMaterialfv(GL_FRONT, GL_EMISSION, material_Ke);
	glMaterialf(GL_FRONT, GL_SHININESS, material_Se);

	// modelview matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//Animation
	teapot();

	// disable lighting
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHTING);

	// swap back and front buffers
	glutSwapBuffers();
}

//================================
// keyboard input
//================================
void keyboard(unsigned char key, int x, int y) {
}

//================================
// reshape : update viewport and projection matrix when the window is resized
//================================
void reshape(int w, int h) {
	// screen size
	g_screenWidth = w;
	g_screenHeight = h;

	// viewport
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);

	// projection matrix
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0, (GLfloat)w / (GLfloat)h, 1.0, 2000.0);
}


//================================
// timer : triggered every 16ms ( about 60 frames per second )
//================================
void timer(int value) {
	// render
	glutPostRedisplay();

	t = t + 0.01;
	if (t >= 1)
	{
		t = 0;
		if (points < number - 4)
		{
			points++;
		}
		else
		{
			points = 0;
		}
	}

	// reset timer
	// 16 ms per frame ( about 60 frames per second )
	glutTimerFunc(16, timer, 0);
}

//================================
// main
//================================
int main(int argc, char** argv) {
	// create opengL window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(1000, 600);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Lab 1");


	// set callback functions
	glutDisplayFunc(render);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutTimerFunc(16, timer, 0);

	// main loop
	glutMainLoop();

	return 0;
}