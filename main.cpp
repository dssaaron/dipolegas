#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>


double part[500][6];
int partAmount = 100;
int mkstepAmount = 200000;
double lambda = 1;
double phi = 0.2;
double xi[100];

double D, d[6];
double Diameter, Length;

double M = 0;

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
	if ((part[N][0]+d[0])*(part[N][0]+d[0])>Length/4) {

		dx = -d[0] + Length - sqrt((2*part[N][0])*(2*part[N][0]));
	}

	//wall collisions
	if ((part[N][1]+d[1])*(part[N][1]+d[1]) + (part[N][2]+d[2])*(part[N][2]+d[2]) > Diameter*Diameter/4) {

		a = d[2]/d[1];
		b = part[N][2] - a*part[N][1];

		if (d[1] > 0) {

			nsX = fmax((-2*a*b - 2*sqrt(Diameter*Diameter/4 * (a*a+1) - b*b))/(2*(a*a+1)), (-2*a*b + 2*sqrt(Diameter*Diameter/4 * (a*a+1) - b*b))/(2*(a*a+1)));
		} else {

			nsX = fmin((-2*a*b - 2*sqrt(Diameter*Diameter/4 * (a*a+1) - b*b))/(2*(a*a+1)), (-2*a*b + 2*sqrt(Diameter*Diameter/4 * (a*a+1) - b*b))/(2*(a*a+1)));
		}

		nsY = a*nsX + b;

		c = part[N][2] + part[N][1]*nsX/nsY;

		dx = 2*(c/(nsY/nsX+nsX/nsY) - part[N][1]);
		dy = 2*(c/(nsY/nsX+nsX/nsY)*nsY/nsX - part[N][2]);
	}


	if ((part[N][1]+dy)*(part[N][1]+dy) + (part[N][2]+dz)*(part[N][2]+dz) < Diameter*Diameter/4) {

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

				distance += (part[N][i]+d[i]-part[partColIterator][i])*(part[N][i]+d[i]-part[partColIterator][i]);
			}

			if (distance <= 1) {

				return false;
			}
		}
	}

	return true;
}

bool energyCheck(int a, double xi) {

	double P[2];
	double distance = 0;
	double e10r, e20r, e11r, e21r, rr1, rr2, e1e20, e1e21;

	//cout << "EC\n";

	P[0] = -xi * part[a][3] * sqrt(1-part[a][5]*part[a][5]);
	P[1] = -xi * d[3] * sqrt(1-d[5]*d[5]);


	for (int b = 0; b < partAmount; b++) {

		if (b != a) {

			for (int i = 0; i < 3; i++) {

				distance += (part[a][i] - part[b][i]) * (part[a][i] - part[b][i]);
			}

			distance = sqrt(distance);

			e10r = (part[a][3] * sqrt(1-part[a][5]*part[a][5]) * (part[b][0]-part[a][0]) +
					part[a][4] * sqrt(1-part[a][5]*part[a][5]) * (part[b][1]-part[a][1]) +
                    part[a][5] * (part[b][2]-part[a][2]));

			e20r = (part[b][3] * sqrt(1-part[b][5]*part[b][5]) * (part[b][0]-part[a][0]) +
					part[b][4] * sqrt(1-part[b][5]*part[b][5]) * (part[b][1]-part[a][1]) +
                    part[b][5] * (part[b][2]-part[a][2]));

			rr1 = distance*distance;

			e1e20 = (part[a][3] * sqrt(1-part[a][5]*part[a][5]) * part[b][3] * sqrt(1-part[b][5]*part[b][5]) +
                     part[a][4] * sqrt(1-part[a][5]*part[a][5]) * part[b][4] * sqrt(1-part[b][5]*part[b][5]) +
                     part[a][5] * part[b][5]);

			P[0] += -lambda * ((3*e10r*e20r/rr1)-e1e20) /(distance*distance*distance);
			
			distance = 0;


			for (int i = 0; i < 3; i++) {
                
            	distance += ((part[a][i]+d[i]) - part[b][i]) * ((part[a][i]+d[i]) - part[b][i]);
        	}

        	distance = sqrt(distance);
            
            e11r = d[3] * sqrt(1-d[5]*d[5]) * (part[b][0] - d[0] - part[a][0]) +
                    d[4] * sqrt(1-d[5]*d[5]) * (part[b][1] - d[1] - part[a][1]) +
                    d[5] * (part[b][2] - d[2] - part[a][2]);
            
            e21r = part[b][3] * sqrt(1-part[b][5]*part[b][5]) * (part[b][0] - d[0] - part[a][0]) +
                    part[b][4] * sqrt(1-part[b][5]*part[b][5]) * (part[b][1] - d[1] - part[a][1]) +
                    part[b][5] * (part[b][2] - d[2] - part[a][2]);
            
            rr2 = distance*distance;
            
            e1e21 = d[3] * sqrt(1-d[5]*d[5]) * part[b][3] * sqrt(1-part[b][5]*part[b][5]) +
                    d[4] * sqrt(1-d[5]*d[5]) * part[b][4] * sqrt(1-part[b][5]*part[b][5]) +
                    d[5] * part[b][5];

            P[1] += -lambda * ((3*e11r*e21r/rr2)-e1e21) /(distance*distance*distance);

		}
	}

	P[0] = exp(P[0]);
    P[1] = exp(P[1]);

	if (P[0] > P[1]) {
        
        return true;
    } else {
        
        if (P[0]/P[1] > random(0, 1)) {
            
            return true;
        } else {
            
            return false;
        }
    }
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
		part[i][0] = random(-Length, Length);
		part[i][1] = random(-Diameter/2, Diameter/2);
		part[i][2] = random(-sqrt(Diameter*Diameter/4 - part[i][1]*part[i][1]), sqrt(Diameter*Diameter/4 - part[i][1]*part[i][1]));

		part[i][3] = cos(angle);
		part[i][4] = sin(angle);
		part[i][5] = random(-1,1);
	}


	//set points where we want to get Magnetizing
	int pointAmount = 1;
	double step;
	for (int i = 0; i < pointAmount; i++) {

		step = i;
		xi[i] = 0.3 + 0*step/pointAmount;
		cout << xi[i] << "; ";
	}
	cout << "\n";

	//directly main cycle
	for (int pointAmountIterator = 0; pointAmountIterator < pointAmount; pointAmountIterator++) {

		for (int mkstepIterator = 0; mkstepIterator < mkstepAmount; mkstepIterator++) {

			if (mkstepIterator%(mkstepAmount/10) == 0) {
            
                cout << "[" << mkstepIterator << "]\n";
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


				if (partCollisionCheck(partAmountInterator) != false) {

					if (energyCheck(partAmountInterator, xi[pointAmountIterator]) != false) {

						for (int i = 0; i < 3; i++) {

							part[partAmountInterator][i] += d[i];
						}

						for (int i = 3; i < 6; i++) {

							part[partAmountInterator][i] = d[i];
						}
					}
				}
			}

			//Magnetizing sum
			if (mkstepIterator > 30000) {

				for (int i = 0; i < partAmount; i++) {

					M += part[i][3]*sqrt(1-part[i][5]*part[i][5]);
				}
			}
		}

		M /= ((mkstepAmount-30000)*M_PI*Diameter*Diameter*Length/4);
		cout << "xi = " << xi[pointAmountIterator] << "; M = " << M << "\n";
		M = 0;
		
	}

	return 0;
}