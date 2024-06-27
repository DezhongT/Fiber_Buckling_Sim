/**
 * simDER
 * simDER stands for "[sim]plified [D]iscrete [E]lastic [R]ods"
 * Dec 2017
 * This code is based on previous iterations.
 * */

//This line is for mac
//#include <GLUT/glut.h>

//This is for linux
#include <GL/glut.h>

#include <iostream>
#include <fstream>
#include <string>
#include "eigenIncludes.h"

// Rod and stepper are included in the world
#include "world.h"
#include "setInput.h"


world myWorld;
int NPTS;
ofstream outfile;
ofstream nodefile;



static void Key(unsigned char key, int x, int y)
{
  switch (key) // ESCAPE to quit
  {
    case 27:
        exit(0);
  }
}

/* Initialize OpenGL Graphics */
void initGL()
{
	glClearColor(0.7f, 0.7f, 0.7f, 0.0f); // Set background color to black and opaque
	glClearDepth(10.0f);                   // Set background depth to farthest
	//glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
	//glDepthFunc(GL_LEQUAL);    // Set the type of depth-test
	glShadeModel(GL_SMOOTH);   // Enable smooth shading
	//glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  // Nice perspective corrections

	glLoadIdentity();
	gluLookAt(0.05, 0.05, 0.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
	glPushMatrix();

	//glMatrixMode(GL_MODELVIEW);
}

void DrawCircle(float cx, float cy, float r, int num_segments)
{
    glBegin(GL_LINE_LOOP);
    for(int ii = 0; ii < num_segments; ii++)
    {
        float theta = 2.0f * 3.1415926f * float(ii) / float(num_segments);//get the current angle

        float x = r * cosf(theta);//calculate the x component
        float y = r * sinf(theta);//calculate the y component

        glVertex3f(x + cx, 0,  y + cy);//output vertex

    }
    glEnd();
}

void DrawCircle2(float cx, float cy, float r, int num_segments)
{
    glBegin(GL_LINE_LOOP);
    for(int ii = 0; ii < num_segments; ii++)
    {
        float theta = 2.0f * 3.1415926f * float(ii) / float(num_segments);//get the current angle

        float x = r * cosf(theta);//calculate the x component
        float y = r * sinf(theta);//calculate the y component

        glVertex3f(x + cx, y + cy, 0);//output vertex

    }
    glEnd();
}

void display(void)
{
	double currentTime  = 0;
  int curr_t = myWorld.currt/1.0 + 1;
	while ( myWorld.simulationRunning() > 0)
	{
		//  Clear screen and Z-buffer
		glClear(GL_COLOR_BUFFER_BIT);

		// draw axis
		double axisLen = 1;
		glLineWidth(0.5);

		glBegin(GL_LINES);
			glColor3f(1.0, 0.0, 0.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(axisLen, 0.0, 0.0);

			glColor3f(0.0, 1.0, 0.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(0.0, axisLen, 0.0);

			glColor3f(0.0, 0.0, 1.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(0.0, 0.0, axisLen);
		glEnd();

    DrawCircle(0, 0, 0.5, 100);
    DrawCircle2(0, 0, 0.5, 100);

		//draw a line
		glColor3f(0.1, 0.1, 0.1);
		glLineWidth(3.0);

    for (int c = 0; c < myWorld.total_rods; c++)
    {
      glBegin(GL_LINES);
      for (int i=0; i < NPTS-1; i++)
      {
        glVertex3f( myWorld.getScaledCoordinate(4*i, c), myWorld.getScaledCoordinate(4*i+1, c), myWorld.getScaledCoordinate(4*i+2, c));
        glVertex3f( myWorld.getScaledCoordinate(4*(i+1), c), myWorld.getScaledCoordinate(4*(i+1)+1, c), myWorld.getScaledCoordinate(4*(i+1)+2, c));
      }
      glEnd();
    }


		glFlush();

    myWorld.updateTimeStep();
    double time = curr_t *1;
    if (abs(time - myWorld.getCurrentTime())< 1e-3)
    {
      myWorld.CoutData(nodefile); // write data to file
      myWorld.CoutDataC(outfile); // write data to file
      curr_t++;
    }
	}

	exit(1);
}

int main(int argc,char *argv[])
{
	setInput inputData;
	inputData = setInput();
	inputData.LoadOptions(argv[1]);
	inputData.LoadOptions(argc,argv);
	//read input parameters from txt file and cmd


	myWorld = world(inputData);
	myWorld.setRodStepper();

	myWorld.OpenFile(outfile, "contact");
  	myWorld.OpenFile(nodefile, "node");


	bool render = myWorld.isRender();
	if (render) // if OpenGL visualization is on
	{
		NPTS = myWorld.numPoints();
		glutInit(&argc,argv);
		glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize (1000, 1000);
		glutInitWindowPosition (500, 500);
		glutCreateWindow ("simDER");
		initGL();
		glutKeyboardFunc(Key);
		glutDisplayFunc(display);
		glutMainLoop();
	}
	else
	{
    int curr_t = myWorld.currt/1.0 + 1;
		while ( myWorld.simulationRunning() > 0)
		{
			myWorld.updateTimeStep(); // update time step

      double time = curr_t *1;
      if (abs(time - myWorld.getCurrentTime())< 1e-3)
      {
        myWorld.CoutData(nodefile); // write data to file
        myWorld.CoutDataC(outfile); // write data to file
        curr_t++;
      }
		}
	}

  // Close (if necessary) the data file
  myWorld.CloseFile(outfile);
  myWorld.CloseFile(nodefile);


	return 0;
}
