#ifndef COLLISION_H
#define COLLISION_H

#include "eigenIncludes.h"
#include "elasticRod.h"
#include "timeStepper.h"

struct edge_pair{
  int node_idx;
  int rod_idx;
  int node_idx1;
  int rod_idx1;
};



struct collisionP{
    Vector3d P; //P is the distance vector between two segments
    double distance;

    double sc; //ratio in segment i -> i+1
    double tc; //ratio in segment j -> j+1

    double kesi;
};

class collision
{
public:
    collision(vector<elasticRod *> m_rod, vector<timeStepper*> m_stepper, int m_rods, double m_radiusofSphere);
    ~collision();

    bool prepareEdges(int curr);
    void computeFc_all();
    void computeFc(int c);
    void setZero();
    double showMinDis();
    int contactNum();
    double getCloestDis();
    int min_i;
    int min_j;
    double computeDis(int i, int j);
    void printInfo(int i, int j);

    void findContactonSphere(int curr);

    void computeWallForce();



private:
    elasticRod *rod;
    timeStepper *stepper;

    vector<elasticRod*> rod_vec;
    vector<timeStepper*> stepper_vec;
    int rods;

    double radiusofSphere;
    double contact_len;

    vector<edge_pair> possibleC;
    vector<edge_pair> contactSphere;
    vector<edge_pair> possibleEnd;


    edge_pair contactP;

    VectorXd Force;

    int ia;

    MatrixXd nodes;

    int nv;

    Vector3d x1s;
    Vector3d x1e;
    Vector3d x2s;
    Vector3d x2e;

    Vector3d d1;
    Vector3d d2;
    Vector3d d12;
    double D1;
    double D2;
    double R;
    double S1;
    double S2;
    double den;
    double t;
    double u;
    double uf;
    Vector3d dist;


    double col_limit;

    collisionP detectionP;

    void computeDistance();

    MatrixXd temp_f;

    vector<double> dis_vec;



};
#endif
