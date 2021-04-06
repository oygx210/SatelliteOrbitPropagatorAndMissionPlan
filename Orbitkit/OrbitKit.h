// OrbitKit.h: interface for the COrbitKit class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ORBITKIT_H_INCLUDED_)
#define AFX_ORBITKIT_H_INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "./Include/Constant.h"


typedef struct	
{
	   double lat;
	   double lon;
	   double alt;
	   double theta;
}geodetic_t;

// General three-dimensional vector structure used by SGP4/SDP4 code. 
typedef struct
{
	double x;
	double y;
	double z;
	double w;
}vector_t;

typedef struct	
{  
	char line1[70];
	char line2[70];
	char name[25];
	long catnum;
	long setnum;
	char designator[10];
	int year;
	double refepoch;
	double incl;
	double raan;
	double eccn;
	double argper;
	double meanan;
	double meanmo;
	double drag;
	double nddot6;
	double bstar;
	long orbitnum;
} sat_t;


typedef struct tag_ORBIT_ELE
{
    double semi_ma;
    double ecc;
	double inc;
    double raan;
    double arg_perigee;
	double omega;		//For SGPSDP prediction model only
    double mean_M0;
	double jd_epoch;
	double jdascnode;	//For SSO
}ORBIT_ELE,*pORBIT_ELE;

typedef struct tag_SAT_ATT
{
	double mass;
	double Ixx;
	double Iyx;
	double Iyy;
	double Izx;
	double Izy;
	double Izz;
	int	rotseq;
	double phi;		//卫星滚动角，度
	double theta;	//卫星俯仰角，度
	double psai;	//卫星偏航角，度
	double q[4];	//欧拉四元素
	double Ispg;	//thruster specific impulse,m/s
	double pulsetm; //thruster's min. pulse time,s
	double appsx;	//X attitude thruster's position w.r.t. CoG,m
	double appsy;	//Y attitude thruster's position w.r.t. CoG,m
	double appsz;	//Z attitude thruster's position w.r.t. CoG,m
	double appsforce;//attitude thruster's force,N
	double oppsx;	//X orbit thruster's position w.r.t. CoG,m
	double oppsy;	//Y orbit thruster's position w.r.t. CoG,m
	double oppsz;	//Z orbit thruster's position w.r.t. CoG,m
	double oppsforce;//orbit thruster's force,N
}SAT_ATT,*pSAT_ATT;

typedef struct tag_SAT_OTHERS
{
	double lifetime;		//the life time of the satellite,year
	double area;			//the surface of the satellite acted by the atmosphere,m2
	double Cd;				//the drag coefficent
}SAT_OTHERS,*pSAT_OTHERS;

typedef struct tag_SAT_POS
{
	double ma;		//卫星半长轴，km
	double ecc;		//卫星偏心率
	double inc;		//卫星轨道倾角，度
	double raan;	//卫星的升交点赤经，度
    double u;		//卫星的纬度幅角，度
	double w;		//卫星的近地点幅角，度
	double geolon;	//地理经度，度
    double geolat;	//地理纬度，度
	double alpha;	//赤经，度
	double delta;	//赤纬，度
	double x;		//三维坐标，km
	double y;		//三维坐标，km
	double z;		//三维坐标，km
	double r;		//半长轴，km
	double vx;		//X方向的速度，km/s
	double vy;		//Y方向的速度，km/s
	double vz;		//Z方向的速度，km/s
}SAT_POS,*pSAT_POS;

//Define the data structure of a single satellite
typedef struct tag_SATELLITE
{
	int id;			// identifier of satellite     
	char modeltype;	// satellite model's type,0=new created,1=SGP4SDP4
	char name[128];	// satellite's name
	bool isactived;	// for SGPSDP model, the satellite's active status
	ORBIT_ELE oe;	// classical orbital elements
	SAT_ATT sa;		// Satellite's current attitude
	SAT_POS sp;		// satellite's current position 
	SAT_OTHERS so;	// satellite's others properties
	sat_t sstle;	// only for SGPSDP model, saving the satellite's orbit parameters
}SATELLITE,*pSATELLITE;

typedef struct tag_XSEGMENT
{
	short 	x1;
	short	y1;
	short	x2;
	short	y2;
}XSegment;

//定义星下点视见区圆的轨迹的点的数量
#define NSEGSVC			360
//定义星下点视见区幅宽线的轨迹的点的数量
#define NSEGSVL			90

const double R_SUN = 696000.0;

class COrbitKit
{
public:
	COrbitKit();
	virtual ~COrbitKit();

//Attributes
public:

//Operations
protected:
	double Kepler(double MeanAnomaly,double Eccentricity);
	double Subtend(double ra1,double dec1,double ra2,double dec2);
	double Zone(short use_dst,double stdz,double jd,double jdb,double jde);
	void Geocent(double geolat,double geolong,double height,double *x_geo,double *y_geo,double *z_geo);
	double Atan22(double x,double y);
	void RotateECL(double jd, double *x, double *y, double *z);
	double CorrectET(double jd);
	void multMatVec(double amVec[3],double pbmVec[3],double mtx[3][3]);
	double absol(double absVec[3]);
	int sign(double data);
	virtual double Latwgs84(double delta = 0);
	virtual double reduce(double value, double rangeMin, double rangeMax);

	void SatOrbitPredict(double jd, SATELLITE *pSat);

public:
	double FixAngle(double angle);
	double FixAngle180(double angle);


	double SSOInc(double alt,double ecc=0);
	double SSOAlt(double inc,double ecc=0);
	double AlphaG( double jd );
	double Raan(double jd);
	void AzimuthElevation(double *Azimuth,double *Elevation,double *Range,double SiteLat,double SiteLong,double SiteAltitude,double SatLat,double SatLong,double SatRadius);
	void MapToOrbitCoord(double *x,double *y,double *z,double alphas,double deltas, double omig,double i,double u);
	void SatellitePosition(double jd, SATELLITE *pSat);
//	void RV2OE(double o_oe[],double i_rv[]);
	void OE2RV(double o_rv[],double i_oe[]);
	double Tsso2Tp(double jdsso,double a,double w=0,double e=0);
	double DeltaRAAN(double dt,double ma,double inc,double ecc);
	double DeltaRAAN(double dt);
	double TNode(double ma,double inc,double ecc=0);
	void SSRO(double *ma,double *inc,double m=30,double n=450,double ecc=0);
	double Adot(double rho,double area,double mass,double a,double Cd=2.5);
	double DeltaA(double rho,double area,double mass,double a,double dL=50,double Cd=2.5);
	double DeltaL(double rho,double area,double mass,double a,double dt=15,double Cd=2.5);

};

#endif
