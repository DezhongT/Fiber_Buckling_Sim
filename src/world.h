#ifndef WORLD_H
#define WORLD_H

#include "eigenIncludes.h"

// include elastic rod class
#include "elasticRod.h"

// include force classes
#include "elasticStretchingForce.h"
#include "elasticBendingForce.h"
#include "elasticTwistingForce.h"
#include "externalGravityForce.h"
#include "inertialForce.h"

// include external force
#include "dampingForce.h"

// include time stepper
#include "timeStepper.h"

// include input file and option
#include "setInput.h"

//include collision checker
#include "collision.h"

extern double* meta_data;


class world
{
public:
	world();
	world(setInput &m_inputData);
	~world();
	void setRodStepper();
	void updateTimeStep();
	int simulationRunning();
	int numPoints();
	double getScaledCoordinate(int i, int c);
	double getCurrentTime();
	double getTotalTime();

	bool isRender();

	// file output
	void OpenFile(ofstream &outfile, string c_name);
	void CloseFile(ofstream &outfile);
	void CoutData(ofstream &outfile);
	void CoutDataC(ofstream &outfile);

	void update(int num);
	double RodLength;
	double radiusofSphere;
	int total_rods;
	double deltaTime;

	double currt;


private:

	// Physical parameters
	double rodRadius;
	int numVertices;
	double youngM;
	double Poisson;
	double shearM;
	double totalTime;
	double density;
	Vector3d gVector;
	double viscosity;

	bool render; // should the OpenGL rendering be included?
	bool saveData; // should data be written to a file?
	double data_period;
	double insertion_time;
	double wait_time;


	// Viscous drag coefficients
	double eta_per, eta_par;

	double tol, stol;
	int maxIter; // maximum number of iterations
	double characteristicForce;
	double forceTol;

	// Geometry
	MatrixXd vertices;
	MatrixXd vertices0;

	VectorXd theta;

	// Rod
	elasticRod *rod;

	vector<elasticRod *> rod_vec;
	vector<timeStepper *> stepper_vec;
	vector<double *> totalForce_vec;



	// set up the time stepper
	timeStepper *stepper;
	double *totalForce;
	double currentTime;

	double time_c;

	// declare the forces
	elasticStretchingForce *m_stretchForce;
	elasticBendingForce *m_bendingForce;
	elasticTwistingForce *m_twistingForce;
	inertialForce *m_inertialForce;
	externalGravityForce *m_gravityForce;
	dampingForce *m_dampingForce;

	//
	vector<elasticStretchingForce *> m_stretchForce_vec;
	vector<elasticBendingForce *> m_bendingForce_vec;
	vector<elasticTwistingForce *> m_twistingForce_vec;
	vector<inertialForce *> m_inertialForce_vec;
	vector<externalGravityForce *> m_gravityForce_vec;
	vector<dampingForce *> m_dampingForce_vec;

	int current_rod;


	collision *check;

	int Nstep;
	int timeStep;
	int iter;

	void rodGeometry(int curr);
	void rodBoundaryCondition();

	void updateBoundary();

	void updateCons();

	void newtonMethod(bool &solved, int iter, bool contact_flag);

	bool calculateForce();

	bool checkcontact();

	double alpha;

	void newtonDamper();

	void updateRod();

	double contact_len;

	double monitor_dis;

	string filename;

};

#endif
