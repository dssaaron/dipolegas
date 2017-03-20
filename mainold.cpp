//
//  main.cpp
//  gastube
//
//  Created by Nick Theone on 12/02/17.
//  Copyright © 2017 Nick Theone. All rights reserved.
//

//
// 1. проверить везде ли общие случаи для размеров частиц и модулей моментов
//

#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <cmath>
#include <cstdlib>

using namespace std;

double M=0;

int N = 100000;                          //mk steps
int partN;                              //part amount

double particle[500][6];                //part array
double partSize = 1;                    //size of particle

double d[6], D;                         //stray & spin

double tubeRatio=10, tubeLen, tubeDiam; //tube

double phi = 0.2;                       //part concentration

double mu = 1.256637E-6;                //
double k = 1.380648E-23;                //
double moment = 1;                      //
double T = 0.7E+16;                         //
double H = 0;                           //

double lambda;
double beta;

double randomPoint(double b, double a) {                                            //selecting random angle cos
    
    return a + rand()%1000000*(b-a)/1000000;
}

double theorMagnet (double H) {
    
    double ksi;
    
    ksi = (moment*mu*H)/(k*T);
    //ksi = mu*moment*((moment * phi *(cosh(ksi)/sinh(ksi) - 1/ksi))/3 + H)/(k*T);
    ksi = H + moment*phi*(cosh(ksi)/sinh(ksi) - 1/ksi)/3;
    return moment * phi *(cosh(ksi)/sinh(ksi) - 1/ksi);
}

bool checkPotential(int a) {
    
    double P[2];
    double distance = 0;
    double e10r, e20r, e11r, e21r, rr1, rr2, e1e20, e1e21;
    double p0, p1;
    
    //part for outside field
    
    //cout << "potential check" << endl;
    
    P[0] = -beta * particle[a][3] * sqrt(1-particle[a][5]*particle[a][5]);
    
    
    P[1] = -beta * d[3] * sqrt(1-d[5]*d[5]);
    
    for (int b = 0; b < partN; b++) {
        
        if (b != a) {
        
        for (int i = 0; i < 3; i++) {
        
            distance += (particle[a][i] - particle[b][i]) * (particle[a][i] - particle[b][i]);
        }
        distance = sqrt(distance);
            
            
            e10r = (particle[a][4] * sqrt(1-particle[a][5]*particle[a][5]) * (particle[b][1]-particle[a][1]) +
                    particle[a][3] * sqrt(1-particle[a][5]*particle[a][5]) * (particle[b][0]-particle[a][0]) +
                    particle[a][5] * (particle[b][2]-particle[a][2]));
            
            e20r = (particle[b][4] * sqrt(1-particle[b][5]*particle[b][5]) * (particle[b][1]-particle[a][1]) +
                    particle[b][3] * sqrt(1-particle[b][5]*particle[b][5]) * (particle[b][0]-particle[a][0]) +
                    particle[b][5] * (particle[b][2]-particle[a][2]));
            
            
            rr1 = distance*distance/partSize/partSize;
            
            e1e20 = (particle[a][3] * sqrt(1-particle[a][5]*particle[a][5]) * particle[b][3] * sqrt(1-particle[b][5]*particle[b][5]) +
                     particle[a][4] * sqrt(1-particle[a][5]*particle[a][5]) * particle[b][4] * sqrt(1-particle[b][5]*particle[b][5]) +
                     particle[a][5] * particle[b][5]);
            
            P[0] += -lambda * ((3*e10r*e20r/rr1)-e1e20) /(distance*distance*distance) * partSize*partSize*partSize;
            
            
        for (int i = 0; i < 3; i++) {
                
            distance += ((particle[a][i]+d[i]) - particle[b][i]) * ((particle[a][i]+d[i]) - particle[b][i]);
        }
        distance = sqrt(distance);
            //cout << "distance new = " << distance << endl;
            
            e11r = d[3] * sqrt(1-d[5]*d[5]) * (particle[b][0] - d[0] - particle[a][0]) +
                    d[4] * sqrt(1-d[5]*d[5]) * (particle[b][1] - d[1] - particle[a][1]) +
                    d[5] * (particle[b][2] - d[2] - particle[a][2]);
            
            e21r = particle[b][3] * sqrt(1-particle[b][5]*particle[b][5]) * (particle[b][0] - d[0] - particle[a][0]) +
                    particle[b][4] * sqrt(1-particle[b][5]*particle[b][5]) * (particle[b][1] - d[1] - particle[a][1]) +
                    particle[b][5] * (particle[b][2] - d[2] - particle[a][2]);
            
            
            rr2 = distance*distance/partSize/partSize;
            
            e1e21 = d[3] * sqrt(1-d[5]*d[5]) * particle[b][3] * sqrt(1-particle[b][5]*particle[b][5]) +
                    d[4] * sqrt(1-d[5]*d[5]) * particle[b][4] * sqrt(1-particle[b][5]*particle[b][5]) +
                    d[5] * particle[b][5];
            
            P[1] += -lambda * ((3*e11r*e21r/rr2)-e1e21) /(distance*distance*distance) * partSize*partSize*partSize;;
            
        
        }
    }
    
    //----------------------
    
    P[0] = exp(P[0]);
    P[1] = exp(P[1]);
    //cout << P[0] << " " << P[1] << endl;
    
    //cout << P[0] << " & " << P[1] << endl;
    
    if (P[0] > P[1]) {
        
        //cout << P[0] << " < " << P[1] << endl;
        return true;
    } else {
        
        if (P[0]/P[1] > randomPoint(1, 0)) {
            
            //cout << P[0]/P[1] << " is lucky\n";
            return true;
        } else {
            
            //cout << P[0]/P[1] << " is unlucky\n";
            return false;
        }
    }
    
    //return true;
}

bool checkJump(int num) {                                                           //check if new position ok
    
    double distance;
    double a, b, c, dx = d[0], dy = d[1], dz = d[2];
    double newSystemx, newSystemy;          //newsystem center

    
    
    // удары об торцы
    
    if ((particle[num][0]+d[0])*(particle[num][0]+d[0]) > tubeLen*tubeLen/4) {
        
        dx = -d[0] + tubeLen - abs(2*particle[num][0]);
        
    }
        
    if (((particle[num][2]+d[2])*(particle[num][2]+d[2])+(particle[num][1]+d[1])*(particle[num][1]+d[1])) > ((tubeDiam)*(tubeDiam)/4)) {    //walls collissions
            
        
        a = d[2]/d[1];
        b = particle[num][2] - a*particle[num][1];
            
        
        if (d[1] > 0) {
            newSystemx = max((-2*a*b - 2*sqrt(pow(tubeDiam/2, 2)*(a*a + 1)-b*b))/(2*(a*a + 1)), (-2*a*b + 2*sqrt(pow(tubeDiam/2, 2)*(a*a + 1)-b*b))/(2*(a*a + 1)));
        } else {
            newSystemx = min((-2*a*b - 2*sqrt(pow(tubeDiam/2, 2)*(a*a + 1)-b*b))/(2*(a*a + 1)), (-2*a*b + 2*sqrt(pow(tubeDiam/2, 2)*(a*a + 1)-b*b))/(2*(a*a + 1)));
        }
            
        
        newSystemy = a*newSystemx+b;
        
            
        c = particle[num][2] + particle[num][1]/newSystemy*newSystemx;
        
            
        dy = 2*(c/(newSystemy/newSystemx+newSystemx/newSystemy) - particle[num][1]);
        dz = 2*(c/(newSystemy/newSystemx+newSystemx/newSystemy)*newSystemy/newSystemx - particle[num][2]);
        
        
    }
    
    if ((particle[num][1] + dy)*(particle[num][1] + dy)+(particle[num][2] + dz)*(particle[num][2] + dz) < (tubeDiam/2)*(tubeDiam/2)) {
        
        d[0] = dx;
        d[1] = dy;
        d[2] = dz;
        
    } else {
        //cout << "doublewall collision \n";
        return false;
    }
    
    
    for (int k = 0; k < partN; k++) {
        if (k != num) {                                                             //another particles collisions
            
            distance = 0;
            
            for (int f = 0; f < 3; f++) {
                
                //distance += pow((particle[num][f]+d[f])-particle[k][f], 2);
                distance +=((particle[num][f]+d[f])-particle[k][f])*((particle[num][f]+d[f])-particle[k][f]);
            }
            //cout << "distance: " << pow(distance, 0.5) << endl;
            if (distance < partSize*partSize) {
                //cout << "particle collision" << endl;
                return false;
            }
        }
    }
    
    
    
    return true;
}

int main() {
    
    int pointN;
    double points[100][2];
    double angle;
    double kappa = 1.948E-1;
    
    srand( static_cast<unsigned int>(time(NULL)));
    
    ofstream output ("data.txt", ios::in | ios::out | ios::app);
    
    cout << "amount of points: \n";
    cin >> pointN;
    
    output << "\n# points: " << pointN << endl;
    
    cout << "points: \n";
    for (int k = 0; k < pointN; k++) {
        
        cin >> points[k][0];
    }
    
    cout << "enter number of particles: ";
    cin >> partN;
    
    output << "# " << partN << " particles" << endl;
    
    lambda = mu*moment*moment/(4*M_PI*partSize*partSize*partSize*k*T);
    
    cout << "\nlambda = " << lambda << endl;
    
    tubeDiam = pow(3*partN*partSize/(2*phi*tubeRatio), 0.3333333333);
    tubeLen = tubeRatio*tubeDiam;
    cout << "\nTube L, D: " << tubeLen << ", " << tubeDiam << "\n" << endl;                   //size of tube by phi coefficient
    
    cout << "kappa = " << kappa << endl;
    
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    for (int i=0; i<partN; i++) {                                                   //generate staring points of particles
        
        angle = randomPoint(2*M_PI, 0);

        particle[i][0] = randomPoint(tubeLen/2, -tubeLen/2);                        // X Y Z
        particle[i][1] = randomPoint(tubeDiam/2, -tubeDiam/2);
        particle[i][2] = randomPoint(sqrt(pow(tubeDiam/2, 2)-pow(particle[i][1],2)), -sqrt(pow(tubeDiam/2, 2)-pow(particle[i][1],2)));
        particle[i][3] = cos(angle);                                        //moment's cos to x
        particle[i][4] = sin(angle);                                        //moment's sin to x
        particle[i][5] = randomPoint(1, -1);                                        //m's cos to z
        
       
    }
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    cout << "\nSUPPOSED RESULTS\n";
    //double tempxi;
    for (int k = 0; k < pointN; k++) {
        
        //tempxi = (moment*mu*points[k][0])/(k*T);
        cout << points[k][0] << " " << theorMagnet(points[k][0]) << endl;
    }
    
    //cout << "\nV = " << M_PI*tubeDiam/4*tubeDiam*tubeLen << endl;
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    for (int pc = 0; pc < pointN; pc++) {
        
        
        H = points[pc][0];
        cout << "\nH = " << H;
        
        beta = (moment*mu*H)/(k*T);
        
        cout << "\nbeta = " << beta << "\n" << endl;
    
        for (int j = 0; j < N; j++) {                                                   //MK steps
        
            if (j%(N/10) == 0) {
            
                cout << "[" << j << "]\n";
            }
        
            for (int i = 0; i < partN; i++) {
            
                angle = randomPoint(2*M_PI, 0);

                D = randomPoint(1, 0);         //stray module
                d[0] = D*randomPoint(1, -1);    //dx
                d[1] = D*randomPoint(1, -1);    //dy ????
                d[2] = D*randomPoint(1, -1);    //dz
                d[3] = cos(angle);               //cos phi for moment
                d[4] = sin(angle);               //sin phi for moment
                d[5] = randomPoint(1, -1);      //cos theta
                
                
            
                if (checkJump(i) != false) {
                
                    if (checkPotential(i) != false) {
                
                        for (int k = 0; k < 3; k++) {       //changing coordinates for particles
                    
                            particle[i][k] += d[k];
                        }
                        
                        for (int k = 3; k < 6; k++) {
                    
                            particle[i][k] = d[k];
                        }
                    
                        
                    } else {
                        
                        //cout << "\nnot jumped\n";
                    }
                    
                    }         //changing particles
            }
        
            if (j > 30000) {
            
                for (int i = 0; i < partN; i++) {
                
                    //M += moment*particle[i][3]*sqrt(1-particle[i][5]*particle[i][5]);
                    points[pc][1] += moment*particle[i][3]*sqrt(1-particle[i][5]*particle[i][5]);
                    //cout << cos(asin(particle[i][3])) << endl;
                }
                //cout << "\n" << M;
            }                         //suming M's
        }
    
        
        points[pc][1] /= ((N-30000)*M_PI*tubeDiam/4*tubeDiam*tubeLen);
        M = 0;
        output << points[pc][0] << " " << points[pc][0] - kappa*points[pc][1] << " " << points[pc][1] << endl;
        
    }       //MK-steps
    
    for (int k = 0; k < pointN; k++) {
        
        cout << "Ho = " << points[k][0] << "; H = " << points[k][0] - kappa*points[k][1] << "; M = " << points[k][1] << endl;
    }
    
    
    
    output.close();
    
    return 0;
}
