#include "collision.h"
#include <bits/stdc++.h>

collision::collision(vector<elasticRod*> m_rod, vector<timeStepper*> m_stepper, int m_rods, double m_radiusofSphere)
{
    rod_vec = m_rod;
    // stepper_vec = m_stepper;
    stepper_vec = m_stepper;
    rods = m_rods;

    radiusofSphere = m_radiusofSphere;

    rod = rod_vec[0];
    Force = VectorXd::Zero(rod->nv * rods * 3);

    ia = ceil(rod->rodRadius * 2 /rod->refLen(0)) + 3;
    nv = rod->nv;

    col_limit = 2;

    contact_len =  2* rod->rodRadius;
    temp_f = MatrixXd::Zero(5,4); //add local force

}

collision::~collision()
{
    ;
}


void collision::findContactonSphere(int curr)
{
  contactSphere.clear();

  for (int c = 0; c < curr; c++)
  {
    rod = rod_vec[c];
    for (int i = 0; i< rod->nv; i++)
    {
      if ((rod->isConstrained[4*i] == 0 && i < rod->nv-1 && rod->isConstrained[4*i+4] == 0)
      || (i == rod->nv - 1 && rod->isConstrained[4*i] == 0))
      {
        Vector3d temp = rod->getVertex(i);

        if (temp.norm() > 0.98 * radiusofSphere)
        {
          contactP.rod_idx = c;
          contactP.node_idx = i;
          contactSphere.push_back(contactP);
        }
      }
    }
  }
}


bool collision::prepareEdges(int curr)
{
  possibleC.clear();
  dis_vec.clear();
  //define the range of nodes to detect collision
  int range = curr * nv;

  nodes = MatrixXd::Zero(range, 3);
  int idx = 0;
  for (int c = 0; c < curr; c++)
  {
    rod = rod_vec[c];
    for (int i = 0; i < rod->nv; i++)
    {
      Vector3d temp = rod->getVertex(i);
      nodes(idx, 0) = temp(0);
      nodes(idx, 1) = temp(1);
      nodes(idx, 2) = temp(2);

      idx++;
    }
  }


  for (int i = 0; i < range - 1 - ia; i++)
  {
    x1s(0) = nodes(i, 0);
    x1s(1) = nodes(i, 1);
    x1s(2) = nodes(i, 2);

    x1e(0) = nodes(i+1, 0);
    x1e(1) = nodes(i+1, 1);
    x1e(2) = nodes(i+1, 2);

    for (int j = i + ia; j < range-1; j++)
    {
      x2s(0) = nodes(j, 0);
      x2s(1) = nodes(j, 1);
      x2s(2) = nodes(j, 2);


      x2e(0) = nodes(j+1, 0);
      x2e(1) = nodes(j+1, 1);
      x2e(2) = nodes(j+1, 2);

      //calculate distance
      computeDistance();
      if (detectionP.distance < contact_len * (1+col_limit))
      {
        contactP.rod_idx = i/nv;
        contactP.node_idx = i%nv;

        contactP.rod_idx1 = j/nv;
        contactP.node_idx1 = j%nv;


        if (rod_vec[contactP.rod_idx]->isConstrained[4*contactP.node_idx]!=2 ||
        rod_vec[contactP.rod_idx1]->isConstrained[4*contactP.node_idx1]!=2 )
        {
          if (contactP.node_idx!=nv-1 && contactP.node_idx1!=nv-1 )
          {
            dis_vec.push_back(detectionP.distance);
            possibleC.push_back(contactP);
            // cout<<contactP.rod_idx<<" "<<contactP.node_idx<<" "<<
            //  contactP.rod_idx1<<" "<<contactP.node_idx1<<" "<<detectionP.distance<<endl;
          }
        }
      }
    }
  }
  // check contact between rod ends;
  for (int i = 0; i < curr-1; i++)
  {
    x1s = rod_vec[i]->getVertex(nv-2);
    x1e = rod_vec[i]->getVertex(nv-1);

    x2s = rod_vec[i+1]->getVertex(0);
    x2e = rod_vec[i+1]->getVertex(1);

    computeDistance();


    if (detectionP.distance < contact_len * (1+col_limit))
    {
      contactP.rod_idx = i;
      contactP.node_idx = nv-2;

      contactP.rod_idx1 = i+1;
      contactP.node_idx1 = 0;


      if (rod_vec[contactP.rod_idx]->isConstrained[4*contactP.node_idx]!=2 ||
      rod_vec[contactP.rod_idx1]->isConstrained[4*contactP.node_idx1]!=2 )
      {
        if (contactP.node_idx!=nv-1 && contactP.node_idx1!=nv-1 )
        {
          dis_vec.push_back(detectionP.distance);
          possibleC.push_back(contactP);
        }
      }
    }
  }

  if (possibleC.size()!=0) return true;

  return false;
}

void collision::computeDistance()
{
  d1 = x1e - x1s;
  d2 = x2e - x2s;
  d12 = x2s - x1s;
  D1 = d1.dot(d1);
  D2 = d2.dot(d2);
  R = d1.dot(d2);
  S1 = d1.dot(d12);
  S2 = d2.dot(d12);

  den = D1*D2 - R*R;

  if (den!= 0)
  {
    t = (S1*D2 - S2 * R)/den;
  }
  else
  {
    t = 0;
  }

  if (t > 1) t = 1;
  if (t < 0) t = 0;
  u = (t * R - S2)/D2;
  uf = u;
  if (uf > 1) uf = 1;
  if (uf < 0) uf = 0;
  if (uf != u) t = (uf * R + S1)/D1;
  if (t > 1) t = 1;
  if (t < 0) t = 0;
  if (uf!= u) u = uf;

  dist = d1 * t - d2 * u - d12;

  detectionP.distance = dist.norm();
  detectionP.P = dist/dist.norm();
  detectionP.sc = t;
  detectionP.tc = u;

}


void collision::setZero()
{
  Force = VectorXd::Zero(rods*nv*3);
  // dis_vec.clear();
}

void collision::computeFc_all()
{
  for (int c =0 ; c < possibleC.size(); c++)
  {
    contactP = possibleC[c];

    x1s = rod_vec[contactP.rod_idx]->getVertex(contactP.node_idx);
    x1e = rod_vec[contactP.rod_idx]->getVertex(contactP.node_idx + 1);

    x2s = rod_vec[contactP.rod_idx1]->getVertex(contactP.node_idx1);
    x2e = rod_vec[contactP.rod_idx1]->getVertex(contactP.node_idx1 + 1);
    computeDistance();
    dis_vec[c] = detectionP.distance;
    if ( contact_len - detectionP.distance > 1e-3 * contact_len)
    {
      detectionP.kesi = (rod_vec[contactP.rod_idx1]->massArray(4*contactP.node_idx1)*(1 - detectionP.tc) +
      rod_vec[contactP.rod_idx1]->massArray(4*(contactP.node_idx1+1))
                              * detectionP.tc)/(rod_vec[contactP.rod_idx]->massArray(4*contactP.node_idx)*(1 - detectionP.sc) +
                               rod_vec[contactP.rod_idx]->massArray(4*(contactP.node_idx+1))*detectionP.sc+ rod_vec[contactP.rod_idx1]->massArray(4*contactP.node_idx)*
                               (1-detectionP.tc) + rod_vec[contactP.rod_idx1]->massArray(4*(contactP.rod_idx1+1))*detectionP.tc);


      temp_f.block<3,1>(0,0) =   detectionP.P * detectionP.kesi * (detectionP.distance - contact_len) * (1-detectionP.sc);
      temp_f.block<3,1>(0,1) =   detectionP.P * detectionP.kesi * (detectionP.distance -  contact_len) * detectionP.sc;
      temp_f.block<3,1>(0,2) =   detectionP.P* (1 - detectionP.kesi) * (contact_len - detectionP.distance) * (1-detectionP.tc);
      temp_f.block<3,1>(0,3) =   detectionP.P* (1 - detectionP.kesi) * (contact_len - detectionP.distance) * detectionP.tc;

      temp_f(3, 0) = contactP.rod_idx;
      temp_f(3, 1) = contactP.rod_idx;
      temp_f(3, 2) = contactP.rod_idx1;
      temp_f(3, 3) = contactP.rod_idx1;



      temp_f(4, 0) = contactP.node_idx;
      temp_f(4, 1) = contactP.node_idx + 1;
      temp_f(4, 2) = contactP.node_idx1;
      temp_f(4, 3) = contactP.node_idx1 + 1;


      for (int c1 = 0; c1 < 4; c1 ++)
      {
        temp_f.block<3,1>(0,c1) = temp_f.block<3,1>(0,c1) *rod_vec[int(temp_f(3, c1))]->massArray(4*int(temp_f(4, c1))) / pow(rod_vec[int(temp_f(3, c1))]->dt, 2);
        for (int c2 = 0; c2 < 3; c2 ++)
        {
           Force(3* (int(temp_f(3,c1)) * nv + int(temp_f(4, c1))) + c2) = Force(3* (int(temp_f(3,c1)) * nv + int(temp_f(4, c1))) + c2) + temp_f(c2, c1);
        }
      }
    }
  }
}

void collision::computeWallForce()
{
  for (int c = 0; c< contactSphere.size(); c++)
  {
    contactP = contactSphere[c];

    Vector3d temp = rod_vec[contactP.rod_idx]->getVertex(contactP.node_idx);

    double d = temp.norm() - radiusofSphere;

    if (d > 0)
    {
      Vector3d force = temp/temp.norm() * d * rod_vec[contactP.rod_idx]->massArray(4*contactP.node_idx)/pow(rod_vec[contactP.rod_idx]->dt,2)*50;

      Force(3*nv * (contactP.rod_idx) + 3*contactP.node_idx) = Force(3*nv * (contactP.rod_idx) + 3*contactP.node_idx) + force(0);
      Force(3*nv * (contactP.rod_idx) + 3*contactP.node_idx + 1) = Force(3*nv * (contactP.rod_idx) + 3*contactP.node_idx + 1)+ force(1);
      Force(3*nv * (contactP.rod_idx) + 3*contactP.node_idx + 2) = Force(3*nv * (contactP.rod_idx) + 3*contactP.node_idx + 2)+ force(2);
    }
  }

}




double collision::showMinDis()
{
  if (dis_vec.size()==0)
  {
    return 1;
  }
  else
  {
    return *min_element(dis_vec.begin(),dis_vec.end());
  }
}

int collision::contactNum()
{
  return dis_vec.size();
}


void collision::computeFc(int c)
{
  int begin = nv * c;
  stepper = stepper_vec[c];

  for (int i = 0; i<nv; i++)
  {
    int j = i + begin;
    Vector3d f_local(Force(3*j), Force(3*j+1), Force(3*j+2));

    for (int k = 0; k<3; k++)
    {
      stepper->addForce(4*i+k, f_local(k));
    }
  }

}


double collision::getCloestDis()
{
  int index = min_element(dis_vec.begin(),dis_vec.end()) - dis_vec.begin();
  contactP = possibleC[index];
  // cout<<index<<endl;
  // cout<<contactP.rod_idx<<" "<<contactP.node_idx<<" "<<
  //  contactP.rod_idx1<<" "<<contactP.node_idx1<<" "<<detectionP.distance<<endl;
  min_i = contactP.rod_idx * nv + contactP.node_idx;
  min_j = contactP.rod_idx1 * nv + contactP.node_idx1;
  return *min_element(dis_vec.begin(),dis_vec.end());
}

double collision::computeDis(int i, int j)
{
  contactP.rod_idx = i/nv;
  contactP.node_idx = i%nv;

  contactP.rod_idx1 = j/nv;
  contactP.node_idx1 = j%nv;

  x1s = rod_vec[contactP.rod_idx]->getVertex(contactP.node_idx);
  x1e = rod_vec[contactP.rod_idx]->getVertex(contactP.node_idx + 1);

  x2s = rod_vec[contactP.rod_idx1]->getVertex(contactP.node_idx1);
  x2e = rod_vec[contactP.rod_idx1]->getVertex(contactP.node_idx1 + 1);

  computeDistance();

  return detectionP.distance;
}

void collision::printInfo(int i, int j)
{
  contactP.rod_idx = i/nv;
  contactP.node_idx = i%nv;

  contactP.rod_idx1 = j/nv;
  contactP.node_idx1 = j%nv;

  x1s = rod_vec[contactP.rod_idx]->getVertex(contactP.node_idx);
  x1e = rod_vec[contactP.rod_idx]->getVertex(contactP.node_idx + 1);

  x2s = rod_vec[contactP.rod_idx1]->getVertex(contactP.node_idx1);
  x2e = rod_vec[contactP.rod_idx1]->getVertex(contactP.node_idx1 + 1);

  cout<<x1s(0) <<" "<< x1s(1) <<" "<<x1s(2) <<" "<<endl;
  cout<<x1e(0) <<" "<< x1e(1) <<" "<<x1e(2) <<" "<<endl;
  cout<<x2s(0) <<" "<< x2s(1) <<" "<<x2s(2) <<" "<<endl;
  cout<<x2e(0) <<" "<< x2e(1) <<" "<<x2e(2) <<" "<<endl;

}
