#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>


double part[500][6];
int partAmount = 50;
int mkstepAmount = 50000;
double lambda = 0;//
double phi = 0.2;
double xi[100];

double D, d[6];
double Diameter, Length;

double M = 0;

double k = 1.3806E-23;
double T = 0.7E+16;
double mu = 1.2566E-6;
double ximult = k*T/mu;

double Mm[500];

double P, p1, p2;

using namespace std;

double random(double a, double b) {

	return a + rand()%10000 * (b-a)/10000;
}


bool partCollisionCheck(int N) {

	double dx = d[0], dy = d[1], dz = d[2];
	double a, b, c;
	double nsX, nsY;

	//cout << "PCC";

	//sidewall collisions
	if ((part[N][2]+d[2])*(part[N][2]+d[2])>Length/4) {

		//dz = -d[2] + Length - sqrt((2*part[N][2])*(2*part[N][2]));
		dz = -dz;
	}

	//wall collisions
	if ((part[N][1]+d[1])*(part[N][1]+d[1]) + (part[N][0]+d[0])*(part[N][0]+d[0]) > Diameter*Diameter/4) {

		a = d[1]/d[0];
		b = part[N][1] - a*part[N][0];

		if (d[0] > 0) {

			nsX = fmax((-2*a*b - 2*sqrt(Diameter*Diameter/4 * (a*a+1) - b*b))/(2*(a*a+1)), (-2*a*b + 2*sqrt(Diameter*Diameter/4 * (a*a+1) - b*b))/(2*(a*a+1)));
		} else {

			nsX = fmin((-2*a*b - 2*sqrt(Diameter*Diameter/4 * (a*a+1) - b*b))/(2*(a*a+1)), (-2*a*b + 2*sqrt(Diameter*Diameter/4 * (a*a+1) - b*b))/(2*(a*a+1)));
		}

		nsY = a*nsX + b;

		c = part[N][1] + part[N][0]*nsX/nsY;

		dx = 2*(c/(nsY/nsX+nsX/nsY) - part[N][0]);
		dy = 2*(c/(nsY/nsX+nsX/nsY)*nsY/nsX - part[N][1]);
	}


	if ((part[N][1]+dy)*(part[N][1]+dy) + (part[N][0]+dx)*(part[N][0]+dx) < Diameter*Diameter/4) {

		d[0] = dx;
		d[1] = dy;
		d[2] = dz;
	} else {

		return false;
	}


	//particle collisions
	double distance;

	for (int partColIterator = 0; partColIterator < partAmount; partColIterator++) {

		if (partColIterator != N) {

			distance = 0;

			for (int i = 0; i < 3; i++) {	

				distance += (part[N][i]+d[i] - part[partColIterator][i]) * (part[N][i]+d[i] - part[partColIterator][i]);
			}

			if (distance <= 1) {

				return false;
			}
		}
	}

	return true;
}


bool energyCheck(int a, double xi) {

	double dist;
	double e10r, e20r, e11r, e21r, rr1, rr2, e1e20, e1e21;

	P = -xi*(part[a][5]-d[5]);

	//change----------------------------------------------------------
	for (int b = 0; b < partAmount; b++) {

		if (b != a) {

			/*for (int i = 0; i < 3; i++) {

				distance += (part[a][i] - part[b][i]) * (part[a][i] - part[b][i]);
			}*/

			dist = sqrt((part[a][0]-part[b][0])*(part[a][0]-part[b][0]) + (part[a][1]-part[b][1])*(part[a][1]-part[b][1]) + (part[a][2]-part[b][2])*(part[a][2]-part[b][2]));

			e10r = (part[a][3] * sqrt(1-part[a][5]*part[a][5]) * (part[b][0]-part[a][0]) +
					part[a][4] * sqrt(1-part[a][5]*part[a][5]) * (part[b][1]-part[a][1]) +
                    part[a][5] * (part[b][2]-part[a][2]));

			e20r = (part[b][3] * sqrt(1-part[b][5]*part[b][5]) * (part[b][0]-part[a][0]) +
					part[b][4] * sqrt(1-part[b][5]*part[b][5]) * (part[b][1]-part[a][1]) +
                    part[b][5] * (part[b][2]-part[a][2]));

			rr1 = dist*dist;

			e1e20 = (part[a][3] * sqrt(1-part[a][5]*part[a][5]) * part[b][3] * sqrt(1-part[b][5]*part[b][5]) +
                     part[a][4] * sqrt(1-part[a][5]*part[a][5]) * part[b][4] * sqrt(1-part[b][5]*part[b][5]) +
                     part[a][5] * part[b][5]);

			// HERE
			P -= lambda * (3*e10r*e20r/rr1 - e1e20) /(dist*dist*dist);
			
			dist = 0;


			/*for (int i = 0; i < 3; i++) {
                
            	distance += ((part[a][i]+d[i]) - part[b][i]) * ((part[a][i]+d[i]) - part[b][i]);
        	}*/

			dist = sqrt((part[a][0]+d[0]-part[b][0])*(part[a][0]+d[0]-part[b][0]) + (part[a][1]+d[1]-part[b][1])*(part[a][1]+d[1]-part[b][1]) + (part[a][2]+d[2]-part[b][2])*(part[a][2]+d[2]-part[b][2])); 
            
            e11r = d[3] * sqrt(1-d[5]*d[5]) * (part[b][0] - d[0] - part[a][0]) +
                    d[4] * sqrt(1-d[5]*d[5]) * (part[b][1] - d[1] - part[a][1]) +
                    d[5] * (part[b][2] - d[2] - part[a][2]);
            
            e21r = part[b][3] * sqrt(1-part[b][5]*part[b][5]) * (part[b][0] - d[0] - part[a][0]) +
                    part[b][4] * sqrt(1-part[b][5]*part[b][5]) * (part[b][1] - d[1] - part[a][1]) +
                    part[b][5] * (part[b][2] - d[2] - part[a][2]);
            
            rr2 = dist*dist;
            
            e1e21 = d[3] * sqrt(1-d[5]*d[5]) * part[b][3] * sqrt(1-part[b][5]*part[b][5]) +
                    d[4] * sqrt(1-d[5]*d[5]) * part[b][4] * sqrt(1-part[b][5]*part[b][5]) +
                    d[5] * part[b][5];

            // AND HERE
            P += lambda * (3*e11r*e21r/rr2 - e1e21)/ (dist*dist*dist);

            dist = 0;
		}
	}
	//change----------------------------------------------------------

    P = exp(P);
    if (P > random(0,1)) {

    	return true;
    } else {

    	return false;
    }
}


double theoreticalMag(double ksi) {
    
    //ksi = ksi*ximult + phi*(cosh(ksi)/sinh(ksi) - 1/ksi)/3;
    return /*phi * */(cosh(ksi)/sinh(ksi) - 1/ksi);
}


int main(int argc, char const *argv[]) {

	srand( static_cast<unsigned int>(time(NULL)));
	

	//calculate Volume
	int tubeRatio = 10;

	Diameter = pow(3*partAmount/(2*phi*tubeRatio), 0.3333);
	//cout << partAmount << " " << phi << " " << tubeRatio << "\n";
	Length = Diameter*tubeRatio;

	cout << "Tube: " << Diameter << " x " << Length << "\n";


	//generate random points: x y z cos sin cos
	double angle;

	for (int i = 0; i < partAmount; i++) {

		angle = random(2*M_PI, 0);
		part[i][0] = random(-Diameter/2, Diameter/2);
		part[i][1] = random(-sqrt(Diameter*Diameter/4 - part[i][0]*part[i][0]), sqrt(Diameter*Diameter/4 - part[i][0]*part[i][0]));
		part[i][2] = random(-Length, Length);

		part[i][3] = cos(angle);
		part[i][4] = sin(angle);
		part[i][5] = random(-1,1);

		cout << "x = " << part[i][0] << " y = " << part[i][1] << " z = " << part[i][2] << " cos phi = " << part[i][3] << " sin phi = " << part[i][4] << " cos theta = " << part[i][5] << endl;
	}


	//set points where we want to get Magnetizing
	int pointAmount = 1;
	double step;
	for (int i = 0; i < pointAmount; i++) {

		step = i;
		xi[i] = 0.3 + 0*step/pointAmount;
		cout << "THEORY:\n";
		cout << "H = " << xi[i]*ximult << "; xi:= " << xi[i] << "; M = " << theoreticalMag(xi[i]) << "\n";
	}
	cout << endl;


	//directly main cycle
	for (int pointAmountIterator = 0; pointAmountIterator < pointAmount; pointAmountIterator++) {

		for (int i = 0; i < partAmount; i++) {

			Mm[i] = 0;
		}

		cout << "[---------------------]\n[";
		for (int mkstepIterator = 0; mkstepIterator < mkstepAmount; mkstepIterator++) {

			if (mkstepIterator%(mkstepAmount/20) == 0) {
            
                //cout << "[" << mkstepIterator << "]\n";
                cout << "+";
                cout.flush();
            } else if (mkstepIterator == mkstepAmount-1) {

            	cout << "+]" << endl;
            }

			for (int partAmountInterator = 0; partAmountInterator < partAmount; partAmountInterator++) {

				//generate way to new statements x y z cos sin cos
				angle = random(2*M_PI, 0);
				D = random(0, 1);
				d[0] = D*random(-1, 1);
				d[1] = D*random(-1, 1);
				d[2] = D*random(-1, 1);

				d[3] = cos(angle);
				d[4] = sin(angle);
				d[5] = random(-1, 1);

				//cout << "dx = " << d[0] << " dy = " << d[1] << " dz = " << d[2] << " cos = " << d[3] << " sin = " << d[4] << " cos = " << d[5] << endl;
				//cout << partAmountInterator << endl;

				if (partCollisionCheck(partAmountInterator) != false) {

					//cout << "x = " << part[partAmountInterator][0] + d[0] << " y = " << part[partAmountInterator][1] + d[1] << " z = " << part[partAmountInterator][2] + d[2] << " cos phi = " << d[3] << " sin phi = " << d[4] << " cos theta = " << d[5] << endl;
					//cout << sqrt((part[partAmountInterator][0] + d[0])*(part[partAmountInterator][0] + d[0]) + (part[partAmountInterator][1] + d[1])*(part[partAmountInterator][1] + d[1])) << " < " << Diameter/2 << endl;
					
					if (energyCheck(partAmountInterator, xi[pointAmountIterator]) != false) {

						//cout << "cos theta = " << d[5] << endl;

						for (int i = 0; i < 3; i++) {

							part[partAmountInterator][i] += d[i];
						}

						for (int i = 3; i < 6; i++) {

							part[partAmountInterator][i] = d[i];
						}
					}
				}

				if (mkstepIterator > 0) {

					Mm[partAmountInterator] += part[partAmountInterator][5]/(mkstepAmount); //-30000
				}

			}
		}

		for (int i = 0; i < partAmount; i++) {

			cout << "Mm " << Mm[i] << endl;
			M += Mm[i];
			//cout << " " << M << endl;
		}

		//M /= (M_PI*Diameter*Diameter*Length/4);
		cout << "xi = " << xi[pointAmountIterator] << "; M = " << M/partAmount << "\n";
	}

	return 0;
}