
#include <iostream>
#include <assert.h>
#include <math.h>

#include <GL/glut.h>

//initialization
int g_screenWidth = 0;
int g_screenHeight = 0;

//number of Points for Spline
static int points = 0; //points index from 0
static int number = 7; //number of points

//time variables Q(t)
static GLfloat t = 0;

//vectors for computing tangent orientation
static GLfloat tangent[3] = { 0 };
static GLfloat binormal[3] = { 0 };
static GLfloat normal[3] = { 0 };
static GLfloat loopIndex = 0;

// The final M matrix for torso, Left Leg and Right Leg
static GLfloat M[16] = { 0 };			//for torso
static GLfloat tempM[3] = { 0 };		//temporate matrix to store the interpolation track (torso position)
static GLfloat Left[16] = { 0 };			//Left Leg
static GLfloat Right[16] = { 0 };		//Right Leg

//Catmul-Rom Spline M Marix
static GLfloat CRSplineM[16] = { -0.5f, 1.0f, -0.5f, 0.0f,	    // Column 1
												   1.5f, -2.5f, 0.0f, 1.0f,       // Column 2
												  -1.5f, 2.0f, 0.5f, 0.0f,        // Column 3
												    0.5f, -0.5f, 0.0f, 0.0f };    // Column 4

//B Spline M Marix
static GLfloat BSplineM[16] = { -1.0f / 6.0f, 3.0f / 6.0f, -3.0f / 6.0f, 1.0f / 6.0f,  // Column 1
												 3.0f / 6.0f, -6.0f / 6.0f, 0.0f / 6.0f, 4.0f / 6.0f,  // Column 2
												-3.0f / 6.0f, 3.0f / 6.0f, 3.0f / 6.0f, 1.0f / 6.0f,   // Column 3
												 1.0f / 6.0f, 0.0f / 6.0f, 0.0f / 6.0f, 0.0f / 6.0f }; // Column 4

//x,y,z in world Cartesian System
static GLfloat points_pos[7][3] = { { 10, 0, -20 },			          //point 1
													{ -10, 0, -20 },				  //point 2
													{ 5, 0, -10 },					  //point 3
													{ 5, 0, -10 },					  //point 4
													{ -5, 0, -5 },					  //point 5
													{ -3, 0, -5 },					  //point 6
													{ 1, 0, -2 } };					  //point 7

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

// Matrix Multiply : Multiply two 4*4 matrix for Leg Animation Function 
void Matrix4Mult4(GLfloat TempLeft[16], GLfloat TempRight[16], GLfloat MResult[16]){
	MResult[0] = TempLeft[0] * TempRight[0] + TempLeft[4] * TempRight[1] + TempLeft[8] * TempRight[2] + TempLeft[12] * TempRight[3];
	MResult[1] = TempLeft[1] * TempRight[0] + TempLeft[5] * TempRight[1] + TempLeft[9] * TempRight[2] + TempLeft[13] * TempRight[3];
	MResult[2] = TempLeft[2] * TempRight[0] + TempLeft[6] * TempRight[1] + TempLeft[10] * TempRight[2] + TempLeft[14] * TempRight[3];
	MResult[3] = TempLeft[3] * TempRight[0] + TempLeft[7] * TempRight[1] + TempLeft[11] * TempRight[2] + TempLeft[15] * TempRight[3];
	MResult[4] = TempLeft[0] * TempRight[4] + TempLeft[4] * TempRight[5] + TempLeft[8] * TempRight[6] + TempLeft[12] * TempRight[7];
	MResult[5] = TempLeft[1] * TempRight[4] + TempLeft[5] * TempRight[5] + TempLeft[9] * TempRight[6] + TempLeft[13] * TempRight[7];
	MResult[6] = TempLeft[2] * TempRight[4] + TempLeft[6] * TempRight[5] + TempLeft[10] * TempRight[6] + TempLeft[14] * TempRight[7];
	MResult[7] = TempLeft[3] * TempRight[4] + TempLeft[7] * TempRight[5] + TempLeft[11] * TempRight[6] + TempLeft[15] * TempRight[7];
	MResult[8] = TempLeft[0] * TempRight[8] + TempLeft[4] * TempRight[9] + TempLeft[8] * TempRight[10] + TempLeft[12] * TempRight[11];
	MResult[9] = TempLeft[1] * TempRight[8] + TempLeft[5] * TempRight[9] + TempLeft[9] * TempRight[10] + TempLeft[13] * TempRight[11];
	MResult[10] = TempLeft[2] * TempRight[8] + TempLeft[6] * TempRight[9] + TempLeft[10] * TempRight[10] + TempLeft[14] * TempRight[11];
	MResult[11] = TempLeft[3] * TempRight[8] + TempLeft[7] * TempRight[9] + TempLeft[11] * TempRight[10] + TempLeft[15] * TempRight[11];
	MResult[12] = TempLeft[0] * TempRight[12] + TempLeft[4] * TempRight[13] + TempLeft[8] * TempRight[14] + TempLeft[12] * TempRight[15];
	MResult[13] = TempLeft[1] * TempRight[12] + TempLeft[5] * TempRight[13] + TempLeft[9] * TempRight[14] + TempLeft[13] * TempRight[15];
	MResult[14] = TempLeft[2] * TempRight[12] + TempLeft[6] * TempRight[13] + TempLeft[10] * TempRight[14] + TempLeft[14] * TempRight[15];
	MResult[15] = TempLeft[3] * TempRight[12] + TempLeft[7] * TempRight[13] + TempLeft[11] * TempRight[14] + TempLeft[15] * TempRight[15];
}

void Normalization(GLfloat N_tempM[3]) {
	GLfloat squar_quaterion = N_tempM[0] * N_tempM[0] + N_tempM[1] * N_tempM[1] + N_tempM[2] * N_tempM[2];
	if (squar_quaterion != 0) // ensure it's not divided by 0
	{
		GLfloat base_quaternion = sqrt(squar_quaterion);
		N_tempM[0] = N_tempM[0] / base_quaternion;
		N_tempM[1] = N_tempM[1] / base_quaternion;
		N_tempM[2] = N_tempM[2] / base_quaternion;
	}
}

// Vector Multiply : two vectors' cross product
void VectorMult(GLfloat TempV1[3], GLfloat TempV2[3], GLfloat VResult[3])
{
	VResult[0] = TempV1[1] * TempV2[2] - TempV1[2] * TempV2[1];
	VResult[1] = TempV1[2] * TempV2[0] - TempV1[0] * TempV2[2];
	VResult[2] = TempV1[0] * TempV2[1] - TempV1[1] * TempV2[0];
}

void Torso_interpolate(GLfloat p_position[7][3], GLfloat SplineM[16])
{
	//T matrix T
	GLfloat TMatrix[4] = { t * t * t, t * t, t, 1 };

	// Set up Tangent T matrix T = {3*t*t,2*t,1,0} 
	GLfloat TangentTMatrix[4] = { 3 * t * t, 2 * t, 1, 0 };

	//Loop to generate the position interpolation track based on 4 points every time
	for (int i = 0; i < 3; i++)
	{
		//changed by timer
		GLfloat GMatrix[4] = { p_position[points][i],
			p_position[(points + 1) % number][i],
			p_position[(points + 2) % number][i],
			p_position[(points + 3) % number][i] };

		tempM[i] = mix(TMatrix, SplineM, GMatrix);
		tangent[i] = mix(TangentTMatrix, SplineM, GMatrix);
	}
	Normalization(tangent);

	if (points == 0 && loopIndex == 0) //loop starts from the beginning
	{
		GLfloat TempVector[3] = { 1, 0, 0 }; //Pick an arbitrary vector V
		Normalization(TempVector);
		VectorMult(tangent, TempVector, normal);
		Normalization(normal);
		VectorMult(normal, tangent, binormal);
		Normalization(binormal);
		loopIndex++;
	}
	else //loop not start from the beginning
	{
		VectorMult(tangent, binormal, normal);
		Normalization(normal);
		VectorMult(normal, tangent, binormal);
		Normalization(binormal);
	}

	// Generate the Interpolation Matrix M
	M[0] = tangent[0];			// column 1 row 1
	M[1] = normal[0];				// .........row 2
	M[2] = binormal[0];			// .........row 3
	M[3] = 0;							// .........row 4

	M[4] = tangent[1];			// column 2 row 1
	M[5] = normal[1];				// .........row 2
	M[6] = binormal[1];			// .........row 3
	M[7] = 0;							// .........row 4

	M[8] = tangent[2];		   // column 3 row 1
	M[9] = normal[2];			   // .........row 2
	M[10] = binormal[2];		   // .........row 3
	M[11] = 0;						   // .........row 4

	M[12] = tempM[0];		  // column 4 row 1
	M[13] = tempM[1];		  // .........row 2
	M[14] = tempM[2];		  // .........row 3
	M[15] = 1;						  // .........row 4
}

void Animation() {
	Torso_interpolate(points_pos, CRSplineM);
	glLoadMatrixf(M);
	glutSolidCube(1.0);
}

void LLeg()
{
	// First Translate Matrix
	GLfloat LT1[16] = { 1, 0, 0, 0,				//column 1
								 0, 1, 0, 0,				//column 2
								 0, 0, 1, 0,				//column 3
								 0, -1, 0, 1 };			//column 4

	//Rotation Matrix by Z axix
	GLfloat LAngle = (sin(4 * 3.14 * t - 3.14 / 2) * 3.14) / 4;     //animate rotation
	GLfloat LT2[16] = { cos(LAngle), sin(LAngle), 0, 0,			   //column 1
						-sin(LAngle), cos(LAngle), 0, 0,					   //column 2
						0, 0, 1, 0,													   //column 3
						0, 0, 0, 1 };													   //column 4 

	// Second Translate Matrix
	GLfloat LT3[16] = { 1, 0, 0, 0,				//column 1
								 0, 1, 0, 0,				//column 2
							 	 0, 0, 1, 0,				//column 3
								 0, 0, 0.3, 1 };			//column 4

	//T3*T2*T1
	Matrix4Mult4(M, LT2, Left);
	Matrix4Mult4(Left, LT1, Left);
	Matrix4Mult4(Left, LT3, Left);

	//Show left leg
	glLoadMatrixf(Left);
	glScalef(0.3f, 2.0f, 0.3f);
	glutSolidCube(1.0);
}

void RLeg()
{
	// First Translate Matrix
	GLfloat RT1[16] = { 1, 0, 0, 0,				//column 1
								  0, 1, 0, 0,				//column 2
								  0, 0, 1, 0,				//column 3
								  0, -1, 0, 1 };			//column 4

	//Rotation Matrix by Z axis
	GLfloat LAngle = (sin(4 * 3.14 * t - 3.14 / 2) * 3.14) / 4;  //animate rotation
	GLfloat RT2[16] = { cos(-LAngle), sin(-LAngle), 0, 0,		//column 1
						-sin(-LAngle), cos(-LAngle), 0, 0,			    //column 2
						0, 0, 1, 0,													//column 3
						0, 0, 0, 1 };													//column 4 

	// Second Translate Matrix
	GLfloat RT3[16] = { 1, 0, 0, 0,				//column 1
								  0, 1, 0, 0,				//column 2
								  0, 0, 1, 0,				//column 3
								  0, 0, -0.3, 1 };		//column 4

	//T3*T2*T1
	Matrix4Mult4(M, RT2, Right);
	Matrix4Mult4(Right, RT1, Right);
	Matrix4Mult4(Right, RT3, Right);

	//Show right leg
	glLoadMatrixf(Right);
	glScalef(0.3f, 2.0f, 0.3f);
	glutSolidCube(1.0);
}

//================================
// timer : triggered every 16ms ( about 60 frames per second )
//================================
void timer(int value) {
	// render
	glutPostRedisplay();

	t = t + 0.005;
	if (t >= 1)
	{
		t = 0;
		if (points < number - 1)
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
	Animation();
	LLeg();
	RLeg();

	// disable lighting
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHTING);

	// swap back and front buffers
	glutSwapBuffers();
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
	gluPerspective(60, (GLfloat)w / (GLfloat)h, 1.0, 20.0);
}

//================================
// keyboard input
//================================
void keyboard(unsigned char key, int x, int y) {
}

//================================
// main
//================================
int main(int argc, char** argv) {
	// create opengL window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(1000, 800);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Lab 2");


	// set callback functions
	glutDisplayFunc(render);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutTimerFunc(16, timer, 0);

	// main loop
	glutMainLoop();

	return 0;
}
