#include "OrbitKit.h"
#include "DateTime.h"
#include <iostream>
#include <fstream>
#include <iomanip> 

using namespace std;

void Demo()
{
	CDateTime m_dateTime;
	COrbitKit m_orbitKit;

	SATELLITE m_sat;
	m_sat.modeltype = NEWCREATEDSAT;
	m_sat.oe.jd_epoch = m_dateTime.GetStartTime(true);
	m_sat.oe.semi_ma = 6846.17879564991;
	m_sat.oe.ecc = 0.0006657;
	m_sat.oe.inc = 97.251;
	m_sat.oe.omega = 139.7753;
	m_sat.oe.arg_perigee = 204.2566;
	m_sat.oe.mean_M0 = 267.5599;
	
	double jd = m_dateTime.GetStartTime(true);
	double start = jd;
	double end = jd + 0.1;
	double step = 1./86400;

	ofstream out("out.txt");

	for (jd = start;jd < end;jd += step)
	{
		m_orbitKit.SatellitePosition(jd, &m_sat);

		out << setprecision(12)<<jd<<" "<< setprecision(9) << m_sat.sp.x << " " << m_sat.sp.y << " " << m_sat.sp.z << " "<<m_sat.sp.vx << " " << m_sat.sp.vy << " " << m_sat.sp.vz << endl;

	}

	out.close();
}



/*SATELLITE ORBITAL DECAY*/
void CalcOrbitLifetime()
{
	double M = 100;//Satellite mass (kg)
	double A = 1;//Satellite area (m^2)
	double H = 300;//Starting height (km)
	double F10 = 70;//Solar Radio Flux (SFU)
	double Ap = 0;//Geomagnetic A index

	cout << " TIME HEIGHT PERIOD MEAN MOTION DECAY" << endl;
	cout << "(days) (km) (mins) (rev/day) (rev/day^2)" << endl;

	//define some values
	double Re = 6378137, Me = 5.98E+24;// 'Earth radius and mass (all SI units)
	double G = 6.67E-11;//'Universal constant of gravitation
	double GMe = G * Me;
	double pi = 3.1415926535897932384626433832795;
	double T = 0, dT = .1;//'time & time increment are in days
	double D9 = dT * 3600 * 24;// 'put time increment into seconds
	double H1 = 10, H2 = H;// 'H2=print height, H1=print height increment
	double R = Re + H * 1000;// 'R is orbital radius in metres
	double P = 2 * pi * sqrt(R * R * R / Me / G);// 'P is period in seconds
	//now iterate satellite orbit with time
	do {
		double SH = (900 + 2.5 * (F10 - 70) + 1.5 * Ap) / (27 - .012 * (H - 200));
		double DN = 6E-10 * exp(-(H - 175) / SH);// 'atmospheric density
		double dP = 3 * pi * A / M * R * DN * D9;// 'decrement in orbital period
		if (H <= H2)// THEN 'test for print
		{
			double Pm = P / 60;
			double MM = 1440 / Pm;
			double nMM = 1440 / ((P - dP)) / 60;// 'print units
			double Decay = dP / dT / P * MM;// 'rev/day/day
			cout << T <<" "<< H << " " << P / 60 << " " << MM << " " << Decay<<endl;// 'do print
			H2 = H2 - H1;// 'decrement print height
		}
		P = P - dP;
		T = T + dT;// 'compute new values
		R = pow((G * Me * P * P / 4 / pi / pi), 1. / 3.);// 'new orbital radius
		H = (R - Re) / 1000;// 'new altitude (semimajor axis)
	} while (H > 180);
	//'now print estimated lifetime of satellite
	cout << "Re-entry after " << T << " days (" << T / 365.25 << " years)" << endl;
}

int main()
{
//	Demo();
	CalcOrbitLifetime();

	return 0;
}