// OrbitKit.cpp: implementation of the COrbitKit class.
//
//////////////////////////////////////////////////////////////////////
#include <cstdio>
#include "OrbitKit.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

COrbitKit::COrbitKit()
{
}

COrbitKit::~COrbitKit()
{

}


int COrbitKit::sign(double data)
{
	return (data>=0?1:-1);
}

double COrbitKit::Zone(short use_dst,double stdz,double jd,double jdb,double jde)
{
	/* Returns Zone time offset when standard time Zone is stdz,
	   when daylight time begins (for the year) on jdb, and ends
	   (for the year) on jde.  This is parochial to the northern
	   hemisphere.  */
	/* Extension -- specifying a negative value of use_dst reverses
	   the logic for the Southern hemisphere; then DST is assumed for
	   the Southern hemisphere summer (which is the end and beginning
	   of the year. */

	if(use_dst == 0) return(stdz);
	else if((jd > jdb) && (jd < jde) && (use_dst > 0)) return(stdz-1.);
	   /* next line .. use_dst < 0 .. for Southern Hemisphere sites. */
	else if(((jd < jdb) || (jd > jde)) && (use_dst < 0)) return(stdz-1.);
	else return(stdz);
}

/* computes the geocentric coordinates from the geodetic
(standard map-type) longitude, latitude, and height.
These are assumed to be in decimal hours, decimal degrees, and
meters respectively.  Notation generally follows 1992 Astr Almanac,
p. K11 */
void COrbitKit::Geocent(double geolat,double geolong,double height,double *x_geo,double *y_geo,double *z_geo)
{

	double denom, C_geo, S_geo;

	geolat = geolat / DEG_IN_RADIAN;
	geolong = geolong / HRS_IN_RADIAN;
	denom = (1. - FLATTEN) * sin(geolat);
	denom = cos(geolat) * cos(geolat) + denom*denom;
	C_geo = 1. / sqrt(denom);
	S_geo = (1. - FLATTEN) * (1. - FLATTEN) * C_geo;
	C_geo = C_geo + height / EQUAT_RAD;  /* deviation from almanac
		       notation -- include height here. */
	S_geo = S_geo + height / EQUAT_RAD;
	*x_geo = C_geo * cos(geolat) * cos(geolong);
	*y_geo = C_geo * cos(geolat) * sin(geolong);
	*z_geo = S_geo * sin(geolat);
}

double COrbitKit::Atan22(double x,double y)
{
	/* returns radian angle 0 to 2pi for coords x, y --
	   get that quadrant right !! */

	double theta;

	if((fabs(x)<1e-8) && (fabs(y)<1e-8)) 
		return(0.);  /* guard ... */

	theta = atan2(y,x);  /* turns out there is such a thing in math.h */
	while(theta < 0.) theta += TWOPI;
	return(theta);
}

/* rotates ecliptic rectangular coords x, y, z to equatorial (all assumed of date.) */
void COrbitKit::RotateECL(double jd, double *x, double *y, double *z)
{
	double incl;
	double /*xpr,*/ypr,zpr;
	double T;

	T = (jd - J2000) / 36525;  /* centuries since J2000 */

	incl = (23.439291 + T * (-0.0130042 - 0.00000016 * T))/DEG_IN_RADIAN;
		/* 1992 Astron Almanac, p. B18, dropping the
		   cubic term, which is 2 milli-arcsec! */
	ypr = cos(incl) * *y - sin(incl) * *z;
	zpr = sin(incl) * *y + cos(incl) * *z;
	*y = ypr;
	*z = zpr;
	/* x remains the same. */
}

double COrbitKit::CorrectET(double jd)
{

	/* Given a julian date in 1900-2100, returns the correction
           delta t which is:
		TDT - UT (after 1983 and before 1998)
		ET - UT (before 1983)
		an extrapolated guess  (after 1998).

	For dates in the past (<= 1998 and after 1900) the value is linearly
        interpolated on 5-year intervals; for dates after the present,
        an extrapolation is used, because the true value of delta t
	cannot be predicted precisely.  Note that TDT is essentially the
	modern version of ephemeris time with a slightly cleaner
	definition.

	Where the algorithm shifts there is an approximately 0.1 second
        discontinuity.  Also, the 5-year linear interpolation scheme can
        lead to errors as large as 0.5 seconds in some cases, though
 	usually rather smaller. */

	double jd1900 = 2415019.5;
	double dates[22];
	double delts[21];  /* can't initialize this look-up table
            with stupid old sun compiler .... */
	double year, delt;
	int i;

	/* this stupid patch for primitive sun C compilers ....
		do not allow automatic initialization of arrays! */

	for(i = 0; i <= 19; i++) dates[i] = 1900 + (double) i * 5.;
	dates[20] = 1998.;  /* the last accurately tabulated one in the
		2000 Almanac ... */

	delts[0] = -2.72;  delts[1] = 3.86; delts[2] = 10.46;
	delts[3] = 17.20;  delts[4] = 21.16; delts[5] = 23.62;
	delts[6] = 24.02;  delts[7] = 23.93; delts[8] = 24.33;
	delts[9] = 26.77;  delts[10] = 29.15; delts[11] = 31.07;
	delts[12] = 33.15;  delts[13] = 35.73; delts[14] = 40.18;
	delts[15] = 45.48;  delts[16] = 50.54; delts[17] = 54.34;
	delts[18] = 56.86;  delts[19] = 60.78; delts[20] = 62.97;

	year = 1900. + (jd - 2415019.5) / 365.25;

	if(year < 1998. && year >= 1900.) {
		i = ((int)year - 1900) / 5;
		delt = delts[i] +
		 ((delts[i+1] - delts[i])/(dates[i+1] - dates[i])) * (year - dates[i]);
	}

	else if (year >= 1998. && year < 2100.)
		delt = 33.15 + (2.164e-3) * (jd - 2436935.4);  /* rough extrapolation */

	else if (year < 1900) {
		printf("CorrectET ... no ephemeris time data for < 1900.\n");
       		delt = 0.;
	}

	else if (year >= 2100.) {
		printf("CorrectET .. very long extrapolation in delta T - inaccurate.\n");
		delt = 180.; /* who knows? */
	}

	return(delt);
}

double COrbitKit::Subtend(double ra1,double dec1,double ra2,double dec2)
{
	/* angle subtended by two directions in the sky. */

	double x1, y1, z1, x2, y2, z2;
	double theta;

	ra1 = ra1 / DEG_IN_RADIAN;
	dec1 = dec1 / DEG_IN_RADIAN;
	ra2 = ra2 / DEG_IN_RADIAN;
	dec2 = dec2 / DEG_IN_RADIAN;
	x1 = cos(ra1)*cos(dec1);
	y1 = sin(ra1)*cos(dec1);
	z1 = sin(dec1);
	x2 = cos(ra2)*cos(dec2);
	y2 = sin(ra2)*cos(dec2);
	z2 = sin(dec2);
   	theta = acos(x1*x2+y1*y2+z1*z2);
	return(theta);
}

/*
功能：根据输入的时间计算卫星在地球惯性坐标系中的赤经和赤纬，以度为单位
*/
void COrbitKit::SatellitePosition(double jd, SATELLITE *pSat)
{
	if(pSat->modeltype ==NEWCREATEDSAT)
	{
		SatOrbitPredict(jd,pSat);
	}
	else 
		return;
}

void COrbitKit::SatOrbitPredict(double jd, SATELLITE *pSat)
{
	
	double jt=jd-pSat->oe.jd_epoch ;

	double ma=pSat->oe.semi_ma ;
	double ecc=pSat->oe.ecc ;
	double inc=RADIAN(pSat->oe.inc) ;	
	double omega=RADIAN(pSat->oe.raan) ;
	double w=RADIAN(pSat->oe.arg_perigee) ;
	double M0=RADIAN(pSat->oe .mean_M0 );

	//calculate the atmosphere density
//	ma=pSat->sp .ma ;
/*	double rho=m_env.AtmosphereDensity (ma -RE);
	double Cd=pSat->so .Cd ;//coefficient of the drag force
	double sqrtau=sqrt(GM*ma)*(1e+6);//convert km/s to m/s
	double AMRatio=pSat->so .area /pSat->sa .mass ;
	double a_dot=-rho*Cd*AMRatio*sqrtau*86400.0/1000.;//km/day
	ma+=a_dot*jt;
*/
	if(fabs(ecc)>0.9999)
	{
		printf("错误!不能计算偏心率大于1的卫星轨道!\r\n");
		return;
	}
	
	double n=sqrt(GM/pow(ma,3))*SEC_IN_DAY;	//求平角速度，弧度/天
	n=n*(1-1.5*J2*pow((RE/ma),2)*(1.5*pow(sin(inc),2)-1)/pow((1-ecc*ecc),1.5));
	double dw = -1.5*n*J2*pow((RE/ma),2)*(2.5*pow(sin(inc),2)-2)/SQR(1-SQR(ecc));
	w+=dw*jt;
	double domega=-1.5*n*J2*pow((RE/ma),2)*cos(inc)/SQR(1-SQR(ecc));
	omega+=domega*jt;
	double M=M0+n*jt;							//求卫星的平近点角
	double f=Kepler(M,ecc);						//求真近点角
	double r=ma*(1-ecc*ecc)/(1+ecc*cos(f));		//求卫星的地心距
	double u=w+f;								//求卫星的纬度幅角

	double sinu=sin(u);
	double cosu=cos(u);
	double sinomega=sin(omega);
	double cosomega=cos(omega);
	double sininc=sin(inc);
	double cosinc=cos(inc);
	double sinw=sin(w);
	double cosw=cos(w);

	double x=r*(cosu*cosomega-sinu*cosinc*sinomega);
	double y=r*(cosu*sinomega+sinu*cosinc*cosomega);
	double z=r*sinu*sininc;
	
	double tmpv=sqrt(GM/(ma*(1-ecc*ecc)));
	double vx=-tmpv*(cosomega*(sinu+ecc*sinw)+sinomega*cosinc*(cosu+ecc*cosw));
	double vy=-tmpv*(sinomega*(sinu+ecc*sinw)-cosomega*cosinc*(cosu+ecc*cosw));
	double vz= tmpv*(sininc*(cosu+ecc*cosw));
//	double PX=cos(omega)*cos(w)-sin(omega)*cos(inc)*sin(w);
//	double PY=sin(omega)*cos(w)+cos(omega)*cos(inc)*sin(w);
//	double PZ=sin(inc)*sin(omega);
//	double QX=-cos(omega)*sin(w)-sin(omega)*cos(inc)*cos(w);
//	double QY=-sin(omega)*sin(w)+cos(omega)*cos(inc)*cos(w);
//	double QZ=sin(inc)*cos(w);
//	double vx=tmpv*(-sin(f)*PX+(cos(f)+ecc)*QX);
//	double vy=tmpv*(-sin(f)*PY+(cos(f)+ecc)*QY);
//	double vz=tmpv*(-sin(f)*PZ+(cos(f)+ecc)*QZ);

	double alpha=0;
	double delta=0;
	if(r>1e-8)
		delta=asin(z/r);
	if(fabs(x)<1e-8) 
	{
		if(y>0)
			alpha=PI/2.0;
		else
			alpha=3.0*PI/2.0;
	}
	else if(x>0)
		alpha=atan(y/x);
	else
		alpha=atan(y/x)+PI;

/*
	double alpha=omega+atan(tan(u)*cos(sat_i));
	double angle=DEGREE(u);
	angle=FixAngle(angle);
	if(angle>=90&&angle<270)
		alpha=PI+alpha;
	double delta=asin(sin(u)*sin(sat_i));
*/
	alpha=DEGREE(alpha);
	delta=DEGREE(delta);

	pSat->sp.ma =ma;
	pSat->sp.inc =FixAngle(DEGREE(inc));
	pSat->sp.raan =FixAngle(DEGREE(omega));
	pSat->sp.u =FixAngle(DEGREE(u));
	pSat->sp.w =FixAngle(DEGREE(w));

	pSat->sp.alpha=FixAngle(alpha);
	pSat->sp.delta=FixAngle(delta);

	double alphaG=AlphaG(jd);
	pSat->sp.geolon =FixAngle180(pSat->sp.alpha -alphaG);
	pSat->sp.geolat =FixAngle180(Latwgs84(pSat->sp.delta)) ;

	pSat->sp.x =x;
	pSat->sp.y =y;
	pSat->sp.z =z;
	pSat->sp.r =r;
	pSat->sp.vx =vx;
	pSat->sp.vy =vy;
	pSat->sp.vz =vz;

}


double COrbitKit::FixAngle (double angle)
{
	double res=0;

	if(angle>=0)
		res=((angle) - 360.0 * (floor((angle) / 360.0)));
	else 
		res=((angle) - 360.0 * (ceil((angle) / 360.0)));

	return res;
}

double COrbitKit::FixAngle180(double angle)
{
	double res=0;

	if(angle>=0)
		res=((angle) - 360.0 * (floor((angle) / 360.0)));
	else 
		res=((angle) - 360.0 * (ceil((angle) / 360.0)));

	if(res>180.0)
		res-=360.0;
	if(res<-180.0)
		res+=360.0;

	return res;
}

/*
功能：根据输入的时间计算格林威治点的赤经，返回值以“度”为单位
*/
double COrbitKit::AlphaG( double jd )
{
  jd -= J2000;      // set relative to 2000.0
  double jdm = jd / TO_CENTURIES;  // convert jd to julian centuries
  double intPart = floor( jd );
  jd -= intPart;
  double rval = 280.46061837 +
                360.98564736629 * jd +
                .98564736629 * intPart +
                jdm * jdm * ( 3.87933e-4 - jdm / 38710000. );

  rval=FixAngle(rval);

  return rval;

/*
	double jd=367*year-7*(year+(month+9)/12)/4+275*month/9+day+1721013.5;
	double djc=(jd-2415020.0)/36525.0;
	double ru=279.69098+36000.76887*djc+0.00039*djc*djc;
	double rt=0.0041780746*(hour*3600+minute*60+second);
	double gst=ru+rt-180;

	gst=FixAngle(gst);
*/
}

/*
	根据输入的时间jd，求对应的太阳同步轨道的RAAN，degree
*/
double COrbitKit::Raan(double jd)
{
	double alphaG=AlphaG(jd);

	double raan=FixAngle(alphaG);

	if(raan<0)
		raan+=360.;

	return raan;

}

/*
	求解开普勒方程
	输入：
	MeanAnomaly--平近点角，弧度
	Eccentricity-偏心率
	输出：
	真近点角，弧度
*/
double COrbitKit::Kepler(double MeanAnomaly,double Eccentricity)
{
	double E;              // Eccentric Anomaly
	double Error;
	double TrueAnomaly;
 
    E = MeanAnomaly;    // Initial guess
    do
        {
        Error = (E - Eccentricity*sin(E) - MeanAnomaly)
                / (1 - Eccentricity*cos(E));
        E -= Error;
        }
   while (fabs(Error) >= 1e-6);
 
    if (fabs(E-PI) < 1e-6)
        TrueAnomaly = PI;
      else
        TrueAnomaly = 2*atan(sqrt((1+Eccentricity)/(1-Eccentricity))
                                *tan(E/2));
    if (TrueAnomaly < 0)
        TrueAnomaly += TWOPI;
 
    return TrueAnomaly;
}

/*
	已知某一天体的赤经和赤纬，计算其在轨道坐标系中的投影
	输入：
	alphas--天体的赤经，弧度
	deltas--天体的赤纬，弧度
	omig----轨道的升交点赤经，弧度
	i-------轨道的倾角，弧度
	u-------轨道的纬度幅角，弧度
	输出：
	x-------天体在轨道坐标系中的向量在X轴中的分量
	y-------天体在轨道坐标系中的向量在Y轴中的分量
	z-------天体在轨道坐标系中的向量在Z轴中的分量
*/
void COrbitKit::MapToOrbitCoord(double *x,double *y,double *z,double alphas,double deltas, double omig,double i,double u)
{
	deltas*=PI/180.0;
	alphas*=PI/180.0;
	omig*=PI/180.0;
	i*=PI/180.0;
	u*=PI/180.0;

	double sindeltas=sin(deltas);
	double cosdeltas=cos(deltas);
	double sinalphas=sin(alphas);
	double cosalphas=cos(alphas);

	double Six=cosdeltas*cosalphas;
	double Siy=cosdeltas*sinalphas;
	double Siz=sindeltas;

	double sinomig=sin(omig);
	double cosomig=cos(omig);
	double sini=sin(i);
	double cosi=cos(i);
	double sinu=sin(u);
	double cosu=cos(u);

	double X1=-sinu*cosomig-cosu*cosi*sinomig;
	double X2=cosu*cosi*cosomig-sinu*sinomig;
	double X3=cosu*sini;

	double Y1=-sini*sinomig;
	double Y2=sini*cosomig;
	double Y3=-cosi;

	double Z1=sinu*cosi*sinomig-cosu*cosomig;
	double Z2=-cosu*sinomig-sinu*cosi*cosomig;
	double Z3=-sinu*sini;

	*x=Six*X1+Siy*X2+Siz*X3;
	*y=Six*Y1+Siy*Y2+Siz*Y3;
	*z=Six*Z1+Siy*Z2+Siz*Z3;
}

void COrbitKit::AzimuthElevation(double *Azimuth,double *Elevation,double *Range,
				double SiteLat,double SiteLong,double SiteAltitude,
				double SatLat,double SatLong,double SatRadius)
{
    double SinSiteLat,sinSiteLong,SinSatLat,SinSatLong;
    double CosSiteLat,cosSiteLong,CosSatLat,CosSatLong;
    double CosBeta,SinBeta;
    double LongDiff,LatDiff;
 
    SiteLong*=PI/180.0;
	SiteLat*=PI/180.0;
	SatLong*=PI/180.0;
	SatLat*=PI/180.0;

	SinSiteLat = sin(SiteLat); sinSiteLong = sin(SiteLong);
    CosSiteLat = cos(SiteLat); cosSiteLong = cos(SiteLong);
    SinSatLat = sin(SatLat); SinSatLong = sin(SatLong);
    CosSatLat = cos(SatLat); CosSatLong = cos(SatLong);
 
    LongDiff = SatLong-SiteLong;
    if (LongDiff < -PI)
        LongDiff += TWOPI;
    if (LongDiff > PI)
        LongDiff -= TWOPI;
 
    //利用球面三角形公式求卫星与观测点之间的地心夹角beta
	CosBeta = SinSiteLat*SinSatLat+CosSiteLat*CosSatLat*cos(LongDiff);
	if(CosBeta<0)
	{
		*Azimuth=0;
		*Elevation=0;
		*Range=0;
		return;
	}
	SinBeta = sqrt(1-SQR(CosBeta));
 
    //利用球面三角形公式求卫星相对观测点的方位角
	*Azimuth = acos((SinSatLat- SinSiteLat*CosBeta)/CosSiteLat/SinBeta);

	if(LongDiff <0)//判断卫星是否在观测点的左边，此时应为大于180度
		*Azimuth=TWOPI-(*Azimuth);

	LatDiff=SatLat-SiteLat;
	if(fabs(LatDiff)<1e-6)//如果在同一条纬度线上，则判断卫星在观测点的左边还是在右边
	{
		if(LongDiff<0)
			*Azimuth=1.5*PI;
		else
			*Azimuth=PI/2.;
	}

	if (fabs(LongDiff)<1e-6)//如果在同一条经度线上，则判断卫星在观测点的北边还是在南边
	{
		if(SiteLat>SatLat)
			*Azimuth =PI;
		else
			*Azimuth =0;
	}

    //转换为以地心为中心
    SiteAltitude += RE*(1-EarthFlat/2+EarthFlat/2*cos(2*SiteLat));

	*Elevation = atan((SatRadius*CosBeta-SiteAltitude)/(SatRadius*SinBeta));
	if(*Elevation<0)//如果仰角小于0，则卫星在地平面以下，此时就将其仰角设为零度
		*Elevation=0;
    //三角形余弦定理
	*Range = sqrt(SQR(SatRadius) + SQR(SiteAltitude)-2*SatRadius*SiteAltitude*CosBeta);

	*Azimuth=(*Azimuth)*180.0/PI;
	*Elevation=(*Elevation)*180.0/PI;
}

/******************************************************************************/
/*                                                                            */
/* multMatVec: multiplies a matrix with a vector in 3D                        */
/*                                                                            */
/******************************************************************************/

void COrbitKit::multMatVec(double amVec[3],double pbmVec[3],double mtx[3][3])
{
    int i;

    for (i = 0; i <= 2; i++)
    {
        pbmVec[i] = mtx[i][0]*amVec[0] + mtx[i][1]*amVec[1] + mtx[i][2]*amVec[2];
    }

}

double COrbitKit::absol(double absVec[3])
{
    double absVal;

    absVal = sqrt(SQR(absVec[0]) + SQR(absVec[1]) + SQR(absVec[2]));

    return(absVal);
}


double COrbitKit::SSOInc(double alt,double ecc)
{
	double a=RE+alt;
	double dlambda=360./365.25*PI/180./86400.0;
	double tmpval=1.0/pow(1-ecc*ecc,2)*pow(RE/a,3.5);
	double coef=-3./2.*J2*sqrt(GM/pow(RE,3))*tmpval ;
	double cosi=dlambda/coef;
	double inc=acos(cosi)*180./PI;

	return inc;
}

double COrbitKit::SSOAlt(double inc,double ecc)
{
	double coef=-9.964;
	double domeg=360./365.25;
	double cosi=cos(inc*PI/180.);
	double RE_a_7_2=domeg*pow(1-ecc*ecc,2)/coef/cosi;
	double RE_a=pow(RE_a_7_2,2./7.);
	double a=RE/RE_a;
	double alt=a-RE;

	return alt;
}

//void COrbitKit::RV2OE(double o_oe[], double i_rv[])
//{
//	float mu=GM;
//	CVector3 rv((float)i_rv[0],(float)i_rv[1],(float)i_rv[2]);
//	CVector3 vv((float)i_rv[3],(float)i_rv[4],(float)i_rv[5]);
//
//// first do some preliminary calculations
//	CVector3 K(0,0,1);
//	CVector3 hv=CrossProduct(rv,vv);
//	CVector3 nv=CrossProduct(K,hv);
//	float n=(float)sqrt(DotProduct(nv,nv));
//	float h2=DotProduct(hv,hv);
//	float v2=DotProduct(vv,vv);
//	float r=(float)sqrt(DotProduct(rv,rv));
//	CVector3 ev=(rv*(v2-mu/r)-vv*DotProduct(rv,vv))/mu;
//	float p=h2/mu;
//
////now compute the oe's
//	float e=(float)sqrt(DotProduct(ev,ev));		//eccentricity
//	float a= p/(1-e*e);							//semimajor axis
//	float i=(float)acos(hv[2]/sqrt(h2));		//inclination
//
//// if e = 0 then omega is undefined and nu is measured from node
//// if i = 0 then Omega is undefined and omega is measured from I
//// if both are zero then  omega and Omega are undefined and nu is measured from I
//
//	double eps=1e-8;
//	double om=0;
//	double Om=0;
//	double nu=0;
//
//	if (e<eps && fabs(i)>eps)     //circular but not equatorial
//	{
//		e=0;
//		om=0;
//		float argacos=nv[0]/n;
//		if (fabs(argacos)>1)
//			argacos=(float)sign(argacos);
//		Om = acos(argacos);		// RAAN
//		if ( nv[1] < 0 )		// fix quadrant
//			Om = 2*PI-Om;
//		
//		argacos=DotProduct(rv,nv)/r/n;
//		if (fabs(argacos)>1)
//			argacos=(float)sign(argacos);
//		nu = acos(argacos);		// true anomaly
//		if ( DotProduct(rv,vv) < 0 )		// fix quadrant
//			nu = 2*PI-nu;
//	}
//
//	if( e>eps && fabs(i)<eps)       //equatorial but not circular
//	{
//		i=0;
//		Om=0;
//		float argacos=ev[0]/e;        // measured from I
//		if (fabs(argacos)>1)
//			argacos=(float)sign(argacos);
//		om = acos(argacos);		// arg of periapsis
//		if ( ev[2] < 0 )		// fix quadrant
//			om = 2*PI-om;
//		argacos = DotProduct(ev,rv)/e/r;
//		if (fabs(argacos)>1)
//			argacos=(float)sign(argacos);
//		
//		nu = acos(argacos);		// true anomaly
//		if ( DotProduct(rv,vv) < 0 )		// fix quadrant
//			nu = 2*PI-nu;
//	}
//	
//	
//	if( e<eps && fabs(i)<eps)     //circular and equatorial
//	{
//		e=0;
//		om=0;
//		i=0;
//		Om=0;
//		float argacos=rv[0]/r;        // measured from I
//		if (fabs(argacos)>1)
//			argacos=(float)sign(argacos);
//		nu = acos(argacos);		// true anomaly
//		if ( vv[0] > 0 )		// fix quadrant
//			nu = 2*PI-nu;
//	}
//
//	if( e>eps && fabs(i)>eps)     //neither circular nor equatorial
//	{
//		float argacos=nv[0]/n;
//		if (fabs(argacos)>1)
//			argacos=(float)sign(argacos);
//
//		Om = acos(argacos);		// RAAN
//		if ( nv[1] < 0 )		// fix quadrant
//			Om = 2*PI-Om;
//		
//		argacos=DotProduct(nv,ev)/n/e;
//		if (fabs(argacos)>1)
//			argacos=(float)sign(argacos);
//
//		om = acos(argacos);		// arg of periapsis
//		if ( ev[2] < 0 )		// fix quadrant
//			om = 2*PI-om;
//		
//		argacos = DotProduct(ev,rv)/e/r;
//		if (fabs(argacos)>1)
//			argacos=(float)sign(argacos);
//			
//		nu = acos(argacos);		// true anomaly
//		if ( DotProduct(rv,vv) < 0 )		// fix quadrant
//			nu = 2*PI-nu;
//	}
//
//	o_oe[0] = a;
//	o_oe[1] = e;
//	o_oe[2] = DEGREE(i);
//	o_oe[3] = DEGREE(Om);
//	o_oe[4] = DEGREE(om);
//	o_oe[5] = DEGREE(nu);
//}

/*
	transform from orbit elements to position and velocity
	Inputs:
	i_oe[0]---ma
	i_oe[1]---ecc
	i_oe[2]---inc
	i_oe[3]---omega
	i_oe[4]---w
	i_oe[5]---u

    Outputs:
	o_rv[0]---x
	o_rv[1]---y
	o_rv[2]---z
	o_rv[3]---vx
	o_rv[4]---vy
	o_rv[5]---vz

*/
void COrbitKit::OE2RV(double o_rv[], double i_oe[])
{
	double ma=i_oe[0];
	double ecc=i_oe[1];
	double inc=RADIAN(i_oe[2]);
	double omega=RADIAN(i_oe[3]);
	double w=RADIAN(i_oe[4]);
	double u=RADIAN(i_oe[5]);
	double f=u-w;
	double r=ma*(1-ecc*ecc)/(1+ecc*cos(f));		//求卫星的地心距
	
	double sinu=sin(u);
	double cosu=cos(u);
	double sinomega=sin(omega);
	double cosomega=cos(omega);
	double sininc=sin(inc);
	double cosinc=cos(inc);
	double sinw=sin(w);
	double cosw=cos(w);

	double x=r*(cosu*cosomega-sinu*cosinc*sinomega);
	double y=r*(cosu*sinomega+sinu*cosinc*cosomega);
	double z=r*sinu*sininc;
	
	double tmpv=sqrt(GM/(ma*(1-ecc*ecc)));
	double vx=-tmpv*(cosomega*(sinu+ecc*sinw)+sinomega*cosinc*(cosu+ecc*cosw));
	double vy=-tmpv*(sinomega*(sinu+ecc*sinw)-cosomega*cosinc*(cosu+ecc*cosw));
	double vz= tmpv*(sininc*(cosu+ecc*cosw));

	o_rv[0]=x;
	o_rv[1]=y;
	o_rv[2]=z;
	o_rv[3]=vx;
	o_rv[4]=vy;
	o_rv[5]=vz;

}

/*
	从升交点地方时转换到过近地点时刻
*/
double COrbitKit::Tsso2Tp(double jdsso,double a,double w,double e)
{
	double jt=jdsso ;				//Epoch时间为过升交点时刻

	w=RADIAN(w) ;

	double n=sqrt(GM/pow(a,3));	//求平角速度，弧度/s
	if(fabs(e)>0.9999)
	{
		printf("错误!不能计算偏心率大于1的卫星轨道!\r\n");
	}
	
	double f=-RADIAN(w);
	double E=2*atan(sqrt((1-e)/(1+e))*tan(f/2));
	double M=E-e*sin(E);

	double tp=jt-M/n/86400.0;

	return tp;
}

/*
	Calculate the delta Omega
	input:
	dt---the period,day
	ma---the semi major axis,km
	inc---the inclination of the orbit,deg
	ecc---eccentric
*/
double COrbitKit::DeltaRAAN(double dt,double ma,double inc,double ecc)
{
	double n=sqrt(GM/pow(ma,3))*SEC_IN_DAY;	//求平角速度，弧度/天
	inc=RADIAN(inc);
	double domega=-1.5*n*J2*pow((RE/ma),2)*cos(inc)/SQR(1-SQR(ecc))*dt;

	domega=DEGREE(domega);

	return domega;
}

/*
	Calculate the delta raan of the Sun sync. orbit
	input:
	dt---deviation of the time in julian date
	ouput:
	the delta raan in degree
*/
double COrbitKit::DeltaRAAN(double dt)
{
	double domega=360/365.25*dt;

	return domega;
}

/*
	Calculate the node period of the orbit
	input:
	ma---semi major axis,km
	inc---inclination of the orbit,degree
	output:
	TN---the node period of the orbit,second
*/
double COrbitKit::TNode(double ma,double inc,double ecc)
{
	inc*=PI/180.;
	double tmp0=1-ecc*ecc;
	double T0=2*PI*sqrt(pow(ma,3)/GM);
	double TN=T0*(1-1.5*J2*pow(RE/(ma*tmp0),2)*(4*cos(inc)*cos(inc)-1));
	return TN;
}

/*
	Calculate the sun synchronous repeating orbit
	input:
	m---the repeat days' number
	n---the repeat cycles' number
	ecc---the ssro's eccentricity
	output:
	ma---the semi major axis,km
	inc---the inclination of the ssro,degree
*/
void COrbitKit::SSRO(double *ma,double *inc,double m,double n,double ecc/*=0*/)
{
	double i0=90.*PI/180.;//initialize the inclination
	double a0=pow(m/n*(sqrt(GM)/OMEGAE),2./3.);//initialize the semi major axis

	double i=i0,a=a0;//first save the initial data
	double di=0,da=0;//used for calculating the deviation

	double tmp0=1-ecc*ecc;
	
	double tol1=1e-5,tol2=1e-6;
	do
	{
		i0=i;
		a0=a;
		double tmp1=m/n-1.5*J2*pow(RE/(a*tmp0),2)*cos(i);
		double tmp2=sqrt(GM)/OMEGAE;
		double tmp3=1-1.5*J2*pow(RE/(a*tmp0),2)*(4*pow(cos(i),2.)-1);
		a=pow(tmp1*tmp2/tmp3,2./3.);
		double cosi=-0.09892*pow(a/RE,7./2.)*pow(tmp0,2.);
		i=acos(cosi);
		di=i-i0;
		da=a-a0;
	}while(fabs(da)>=tol1&&fabs(di)>=tol2);

	*ma=a;
	*inc=i*180./PI;
}

/*
	Calculate the semi major axis decaying due to the atmosphere drag
	rho---atmosphere density,kg/m3
	area---area of the satellite surface,m2
	mass---mass of the satellite,kg
	a---semi major axis of the satellite,km
	Cd---the drag coefficent
	output:
	the decay rate of the satellite' semi major axis,km/s
*/
double COrbitKit::Adot(double rho,double area,double mass,double a,double Cd/*=2.5*/)
{
	double adot=-rho*Cd*area/mass*sqrt(a*GM*(1e+12))/1000;//km/s

	return adot;

}

/*
	Calculate the delta axis due to the maintenance of the ground track
	rho---atmosphere density,kg/m3
	area---area of the satellite surface,m2
	mass---mass of the satellite,kg
	a---semi major axis of the satellite,km
	dL---the required distance between ground tracks to be adjusted,km
	Cd---the drag coefficent
	output:
	the delta axis to adjusted,km
*/
double COrbitKit::DeltaA(double rho,double area,double mass,double a,double dL/*=50*/,double Cd/*=2.5*/)
{

	double adot= Adot(rho,area,mass,a,Cd)*86400;//km/day
	double deltaA=sqrt(-2*a*adot*dL/(3*PI*RE));
	
	return deltaA;
}

/*
	Calculate the delta longitude due to the decaying of the orbit
	rho---atmosphere density,kg/m3
	area---area of the satellite surface,m2
	mass---mass of the satellite,kg
	a---semi major axis of the satellite,km
	dt---the specified time to be calculated,day
	Cd---the drag coefficent
	output:
	the delta longitude that may be introduced by the decaying of the orbit,km
*/
double COrbitKit::DeltaL(double rho,double area,double mass,double a,double dt/*=15*/,double Cd/*=2.5*/)
{
	double adot= Adot(rho,area,mass,a,Cd)*86400;//km/day
	double deltaL=3*PI*RE*adot*pow(dt,2)/(2*a);

	return deltaL;
}

/*
	计算WGS84坐标系下的地理纬度
	输入：lat---地心纬度值，度
	输出：地理纬度值，度
*/
double COrbitKit::Latwgs84(double delta)
{
	const double f = 1. / 298.257;
	delta = RADIAN(delta);

	double tandelta = tan(delta);
	double tanlat = tandelta / pow(1 - f, 2);
	double geolat = atan(tanlat);

	return DEGREE(geolat);
}

/*
	To limit the value in the range between rangeMin and rangeMax
*/
double COrbitKit::reduce(double value, double rangeMin, double rangeMax)
{
	double range, rangeFrac, fullRanges, retval;

	range = rangeMax - rangeMin;
	rangeFrac = (rangeMax - value) / range;

	modf(rangeFrac, &fullRanges);

	retval = value + fullRanges * range;

	if (retval > rangeMax)
		retval -= range;

	return(retval);
}