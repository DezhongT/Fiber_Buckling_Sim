#include "world.h"
#include <sstream>
#include <iomanip>
#include <float.h>


bool IsFiniteNumber(double x)
{
	return (x <= DBL_MAX && x >= -DBL_MAX);
}

world::world()
{
	;
}

world::world(setInput &m_inputData)
{
	render = m_inputData.GetBoolOpt("render");				// boolean
	saveData = m_inputData.GetBoolOpt("saveData");			// boolean

	// Physical parameters
	RodLength = m_inputData.GetScalarOpt("RodLength");      // meter
  	gVector = m_inputData.GetVecOpt("gVector");             // m/s^2
 	maxIter = m_inputData.GetIntOpt("maxIter");             // maximum number of iterations
	rodRadius = m_inputData.GetScalarOpt("rodRadius");      // meter
	numVertices = m_inputData.GetIntOpt("numVertices");     // int_num
	youngM = m_inputData.GetScalarOpt("youngM");            // Pa
	Poisson = m_inputData.GetScalarOpt("Poisson");          // dimensionless
	deltaTime = m_inputData.GetScalarOpt("deltaTime");      // seconds
	totalTime= m_inputData.GetScalarOpt("totalTime");       // seconds
	tol = m_inputData.GetScalarOpt("tol");                  // small number like 10e-7
	stol = m_inputData.GetScalarOpt("stol");				// small number, e.g. 0.1%
	density = m_inputData.GetScalarOpt("density");          // kg/m^3
	viscosity = m_inputData.GetScalarOpt("viscosity");      // viscosity in Pa-s

	data_period = m_inputData.GetScalarOpt("data-period");
	insertion_time = m_inputData.GetScalarOpt("insertion-time");
	wait_time = m_inputData.GetScalarOpt("wait-time");
	total_rods = m_inputData.GetIntOpt("total-rods");
	boundary_iter = m_inputData.GetIntOpt("boundary-iter");
	external_flow = m_inputData.GetScalarOpt("external-flow");
	flow_period = m_inputData.GetScalarOpt("flow-period");
	radiusofSphere = m_inputData.GetScalarOpt("radiusofSphere");

	filename = m_inputData.GetStringOpt("filename");
	

	shearM = youngM/(2.0*(1.0+Poisson));					// shear modulus

	// Viscous drag coefficients using Resistive Force Theory
	eta_per = 4.0*M_PI*viscosity/( log(2*RodLength/rodRadius) + 0.5);
  	eta_par = 2.0*M_PI*viscosity/( log(2*RodLength/rodRadius) - 0.5);

	contact_len = 2 * rodRadius;

  	alpha = 1;
}

world::~world()
{
	;
}

bool world::isRender()
{
	return render;
}

void world::OpenFile(ofstream &outfile, string c_name)
{
	if (saveData==false) return;

	int systemRet = system("mkdir datafiles"); //make the directory
	if(systemRet == -1)
	{
		cout << "Error in creating directory\n";
	}

	time_t current_time = time(0);

	// Open an input file named after the current time
	ostringstream name;

  name << "datafiles/simDER_" <<to_string(youngM)<<"_"<<to_string(RodLength)<<"_"
	<<to_string(rodRadius)<<"_"<<to_string(total_rods)<<"_"<< c_name<<".txt";
	outfile.open(name.str().c_str());
	outfile.precision(10);
}


void world::CloseFile(ofstream &outfile)
{
	if (saveData==false)
		return;

	outfile.close();
}

void world::CoutData(ofstream &outfile)
{
	if (saveData==false)
		return;


	for (int c = 0; c < total_rods; c++)
	{
		rod = rod_vec[c];
		for (int i =0 ; i<numVertices; i++)
		{
			if (i < numVertices -1)
			{
				outfile<<rod->x(4*i)<<" "<<rod->x(4*i+1)<<" "<<rod->x(4*i+2)<<" "<<rod->x(4*i+3)<<
				" "<<rod->u(4*i)<<" "<<rod->u(4*i+1)<<" "<<rod->u(4*i+2)<<" "<<rod->u(4*i+3)<<
				" "<<rod->tangent(i, 0)<<" "<<rod->tangent(i, 1)<<" "<<rod->tangent(i, 2)<<
				" "<<rod->refTwist(i)<<" "<<rod->d1(i, 0)<<" "<<rod->d1(i, 1)<<" "<<rod->d1(i,2)
				<<" "<<rod->d2(i, 0)<<" "<<rod->d2(i,1)<<" "<<rod->d2(i,2)<<endl;
			}
			else
			{
				outfile<<rod->x(4*i)<<" "<<rod->x(4*i+1)<<" "<<rod->x(4*i+2)<<" "<<0<<
				" "<<rod->u(4*i)<<" "<<rod->u(4*i+1)<<" "<<rod->u(4*i+2)<<" "<<0<<
				" "<<0<<" "<<0<<" "<<0<<
				" "<<0<<" "<<0<<" "<<0<<" "<<0
				<<" "<<0<<" "<<0<<" "<<0<<endl;
			}
		}
	}
}

void world::CoutDataC(ofstream &outfile)
{
	if (saveData==false)
		return;

	outfile<< currentTime<<" "<<check->contactNum()<<endl;
}


void world::setRodStepper()
{
	//total rod
	rod_vec.clear();
	stepper_vec.clear();
	for (int i = 0; i< total_rods; i++)
	{
		// elasticRod *temp;
		rodGeometry(i);
		rod = new elasticRod(vertices, vertices0, density, rodRadius, deltaTime,
			youngM, shearM, RodLength, theta);
	  	rod_vec.push_back(rod);
	}

	// Find out the tolerance, e.g. how small is enough?
	characteristicForce = M_PI * pow(rodRadius ,4)/4.0 * youngM / pow(RodLength, 2);
	forceTol = tol * characteristicForce;

	// Set up boundary condition
	rodBoundaryCondition();

	// setup the rod so that all the relevant variables are populated
	for (int c = 0; c <  total_rods; c++)
	{
		rod = rod_vec[c];
		rod->setup();
	}
	// End of rod setup

	// set up the time stepper
	for (int i = 0; i< total_rods; i++)
	{
		stepper = new timeStepper(* rod_vec[i]);
		stepper_vec.push_back(stepper);
	}

	for (int c = 0; c <  total_rods; c++)
	{
		stepper = stepper_vec[c];
		totalForce = stepper->getForce();
		totalForce_vec.push_back(totalForce);
	}

	for (int c = 0; c <  total_rods; c++)
	{
		stepper = stepper_vec[c];
    	rod = rod_vec[c];
		m_stretchForce = new elasticStretchingForce(*rod, *stepper);
		m_bendingForce = new elasticBendingForce(*rod, *stepper);
		m_twistingForce = new elasticTwistingForce(*rod, *stepper);
		m_inertialForce = new inertialForce(*rod, *stepper);
		m_gravityForce = new externalGravityForce(*rod, *stepper, gVector);
		m_dampingForce = new dampingForce(*rod, *stepper, viscosity, eta_per, eta_par);

		m_stretchForce_vec.push_back(m_stretchForce);
   		m_bendingForce_vec.push_back(m_bendingForce);
		m_twistingForce_vec.push_back(m_twistingForce);
		m_inertialForce_vec.push_back(m_inertialForce);
		m_gravityForce_vec.push_back(m_gravityForce);
		m_dampingForce_vec.push_back(m_dampingForce);
	}
	
	check = new collision(rod_vec, stepper_vec, total_rods, radiusofSphere);

	Nstep = totalTime/deltaTime;

	// Allocate every thing to prepare for the first iteration
	current_rod = 0;

	for (int c = 0; c <  total_rods; c++)
	{
		rod = rod_vec[c];
		rod->updateTimeStep();
	}

	timeStep = 0;
	currentTime = 0;

	cout <<" Sim environment is set up"<<endl;
}

// Setup geometry
void world::rodGeometry(int curr)
{
	 vertices = MatrixXd(numVertices, 3);
	 vertices0 = MatrixXd(numVertices, 3);

	 double offset = 0.01;
	 if (curr == 0) offset = 0;

	 double delta_l = RodLength/(numVertices-1);
	 for (int i = 0; i< numVertices; i++)
	 {
		 vertices(i, 0) = i* delta_l + radiusofSphere + offset;
		 vertices(i, 1) = 1.0*rand()/RAND_MAX * delta_l * 1e-2 ;
		 vertices(i, 2) = 1.0*rand()/RAND_MAX * delta_l * 1e-2;

		 vertices0(i, 0) = i* delta_l + radiusofSphere;
		 vertices0(i, 1) = 0 ;
		 vertices0(i, 2) = 0;
	 }

	 theta = VectorXd::Zero(numVertices - 1);
}


void world::rodBoundaryCondition()
{
	// set up boundary
	for (int c = 0; c< total_rods; c++)
	{
		rod = rod_vec[c];

		for (int i = 0; i< numVertices; i++)
		{
			rod->setVertexBoundaryCondition(rod->getVertex(i), i, 2);
			if (i < numVertices-1) rod->setThetaBoundaryCondition(rod->getTheta(i), i, 2);
		}
	}
}


void world::update(int num)
{

	rod = rod_vec[num];
	stepper = stepper_vec[num];
	totalForce = totalForce_vec[num];
	rod->updateMap();
	stepper->update();
	totalForce = stepper->getForce();
	// rod_vec[num] = rod;
	// stepper_vec[num] = stepper;
	totalForce_vec[num] = totalForce;
}

void world::updateBoundary()
{
	// check if we need to inject the next rod
	rod = rod_vec[current_rod];
	if (rod->isConstrained[4*(numVertices-1)] == 0 && current_rod < total_rods-1)
	{
		current_rod = current_rod + 1;
	}
	rod = rod_vec[current_rod];

	// inject the rod if there exsit parts of the rod outside the tank
	if (rod->isConstrained[4*(numVertices-1)]==2)
	{
	    double u  = RodLength/insertion_time;
		for (int i = 0; i< numVertices; i++)
		{
			//change velocity
			if (rod->x(4*i) >= radiusofSphere)
			{
				rod->setVelocity(i, Vector3d(-u,0,0));
			}
		}
	}
}


void world::updateRod()
{
	rod = rod_vec[current_rod];

	for (int i = 0; i< numVertices; i++)
	{
		Vector3d temp = rod->getVertex(i);
		if (temp(0) < radiusofSphere)
		{
			rod->setVertexBoundaryCondition(temp, i, 0);
			if (i > 0) rod->setThetaBoundaryCondition(0, i-1, 0);
		}
		else
		{
			if (rod->isConstrained[4*i]==0)
			{
				if (temp(1) == vertices(i, 1))
				{
					rod->setVertexBoundaryCondition(temp, i, 2);
				}
			}
		}
	}
	update(current_rod);
}

void world::updateTimeStep()
{
	// update Boundary (in the form of updating velocity
	updateBoundary();

	//step 1: check if we need to use smaller timestep
	while (true)
	{
		//update Guess
		for (int i = 0; i < current_rod + 1; i++)
		{
			rod = rod_vec[i];
			rod->updateGuess();
		}
    	updateRod();

		//sub step 1: check the contact exits in rod segment and rod ends
		bool contact_flag = check->prepareEdges(current_rod+1);
		//sub step 2: check if contact with the tank
		check->findContactonSphere(current_rod + 1);

		// we need to find the guess minimum distance to determine shrink time step size or not
		bool change_time_step = false;
		if (contact_flag)
		{
			double min_dis = check->getCloestDis();
			// cout<<"min i: "<<check->min_i<<" min_j: "<<check->min_j<< " guess min dis: "<<min_dis<<endl;
			if (min_dis < 0.8 * contact_len)
			{
				change_time_step = true;
				deltaTime = 0.1*deltaTime;
				for (int i = 0; i< total_rods; i++)
				{
					rod = rod_vec[i];
					m_dampingForce = m_dampingForce_vec[i];
					rod->dt = deltaTime;
					m_dampingForce->dt = deltaTime;
				}
			}
			else
			{
				break;
			}
		}
		else
		{
			break;
		}
	}

	//Step 2 solving the nonlinear system
	bool solved = false;
	int iter_c = 0;
	bool contact_flag = false;
	while (!solved)
	{
		for (int i = 0; i < current_rod + 1; i++)
		{
			rod = rod_vec[i];
			rod->updateGuess();
		}

		updateRod();
		contact_flag = check->prepareEdges(current_rod+1);
		newtonMethod(solved, 0, contact_flag);

		if (!solved)
		{
			deltaTime = 0.1*deltaTime;
			for (int i = 0; i< total_rods; i++)
			{
				rod = rod_vec[i];
				m_dampingForce = m_dampingForce_vec[i];
				rod->dt = deltaTime;
				m_dampingForce->dt = deltaTime;
			}
		}
		iter_c++;

		if (iter_c>100 || deltaTime < 1e-9)
		{
			currentTime = totalTime;
			break;
		}
	}

	for (int i = 0; i < current_rod + 1; i++)
	{
		rod = rod_vec[i];
		rod->updateTimeStep();
	}

	if (render) cout << "time: " << currentTime << " iter= " << iter <<" deltaTime: "<<deltaTime<<
	" current rod: "<<current_rod<<endl;
	currentTime += deltaTime;


	//Step 4: whether to update time step size or not
	if (iter <= 5 && deltaTime < 5e-3)
	{
		deltaTime = 10*deltaTime;
		for (int i = 0; i< total_rods; i++)
		{
			rod = rod_vec[i];
			m_dampingForce = m_dampingForce_vec[i];
			rod->dt = deltaTime;
			m_dampingForce->dt = deltaTime;
		}
	}

	if (iter >= 20 && deltaTime > 5e-3)
	{
		deltaTime = 0.1*deltaTime;
		for (int i = 0; i< total_rods; i++)
		{
			rod = rod_vec[i];
			m_dampingForce = m_dampingForce_vec[i];
			rod->dt = deltaTime;
			m_dampingForce->dt = deltaTime;
		}
	}

 	timeStep++;

	if (solved == false)
	{
		currentTime = totalTime; // we are exiting
	}
}


void world::newtonDamper()
{
    if (iter < 2)
        alpha = 1.0;
    else if (iter < 4)
        alpha = 0.5;
    else if (iter < 6)
        alpha = 0.25;
    else if (iter < 8)
        alpha = 0.10;
    else if (iter < 10)
        alpha = 0.05;
    else if (iter > 50)
        alpha = 0.01;
}




void world::newtonMethod(bool &solved, int m_iter, bool contact_flag)
{
	double normf = forceTol * 10.0;
	double normf0 = 0;
	iter = m_iter;

	while (solved == false )
	{

		if (iter == 0)
		{
			check->setZero();
		}

		if (contact_flag)
		{
			check->computeFc_all();
		}

		check->computeWallForce();

    	//calculate Force
		normf = 0;
		for (int c = 0; c< current_rod +1; c++)
		{

			rod = rod_vec[c];
			if (rod->isConstrained[0]!= 0) continue;

			stepper = stepper_vec[c];
			m_inertialForce = m_inertialForce_vec[c];
			m_stretchForce = m_stretchForce_vec[c];
			m_bendingForce = m_bendingForce_vec[c];
			m_twistingForce = m_twistingForce_vec[c];
			m_gravityForce = m_gravityForce_vec[c];
			m_dampingForce = m_dampingForce_vec[c];
			totalForce = totalForce_vec[c];

			rod->prepareForIteration();

			stepper->setZero();

		  	// Compute the forces and the jacobians
			m_inertialForce->computeFi();
			m_inertialForce->computeJi();

		  	m_stretchForce->computeFs();
		  	m_stretchForce->computeJs();

		 	m_bendingForce->computeFb();
		  	m_bendingForce->computeJb();

		  	m_twistingForce->computeFt();
		  	m_twistingForce->computeJt();

		  	m_gravityForce->computeFg();
		  	m_gravityForce->computeJg();

		 	m_dampingForce->computeFd();
		  	m_dampingForce->computeJd();

			check->computeFc(c);


		  	// Compute norm of the force equations.
			for (int i=0; i < rod->uncons; i++)
			{

				normf += totalForce[i] * totalForce[i];
			}
		}
		normf = sqrt(normf);

	  	bool checker = false;

		if (contact_flag)
		{
			double min = check->getCloestDis();
			if (min > 0.99 * contact_len)
			{
				checker = true;
			}
			monitor_dis = min;
		}
		else
		{
			checker = true;
		}

		if (iter == 0)
		{
			normf0 = normf;
		}
		else if(iter > 0 && normf <= forceTol *sqrt(current_rod+1) && checker)
		{
			if (IsFiniteNumber(normf)) solved = true;
		}
		else{
			for (int c = 0; c< current_rod +1; c++)
			{

				rod = rod_vec[c];
				if (rod->isConstrained[0]!= 0) continue;

				stepper = stepper_vec[c];
				totalForce = totalForce_vec[c];
			

				if (solved == false)
				{
					stepper->integrator(); // Solve equations of motiona = 1;
					rod->updateNewtonX(stepper->DX, alpha); // new q = old q + Delta q
				}
			}
		}
		if (iter > maxIter)
		{
			if (normf < 1e-10)
			{
				solved = true;
			}
			cout << "Error. Could not converge. adjusting time step size .\n";
			break;
		}
		iter++;
	}
}

int world::simulationRunning()
{
	if (currentTime<totalTime)
		return 1;
	else
	{
		return -1;
	}
}

int world::numPoints()
{
	return rod->nv;
}

double world::getScaledCoordinate(int i, int c)
{
	rod = rod_vec[c];
	if (rod->isConstrained[i] != 2)
	{
		return rod->x[i] /(radiusofSphere * 2);
	}
	else
	{
		if (i%4 == 0)
		{
			return radiusofSphere/(radiusofSphere * 2);
		}
		return 0;
	}
}

double world::getCurrentTime()
{
	return currentTime;
}

double world::getTotalTime()
{
	return totalTime;
}
