// DateTime.cpp: implementation of the CDateTime class.
//
//////////////////////////////////////////////////////////////////////
#include <windows.h>
#include "DateTime.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CDateTime::CDateTime()
{
	m_year=2021;
	m_month=1;
	m_day=1;
	m_hour=0;
	m_minute=0;
	m_second=0;

	m_epochjd=JD(m_year,m_month,m_day,m_hour,m_minute,m_second);
	m_startjd= m_epochjd;
	m_currentjd= m_epochjd;
	m_endjd= m_epochjd +0.5;
	m_stepjd=1.0/86400.0;

}

CDateTime::~CDateTime()
{

}

/*
	根据输入的时间计算Julian时间，返回值以“天”为单位
*/
double CDateTime::JD(long year, int month, int day, int hour, int minute, int second)
{
	double jd;
 	int y, m;
	long A, B;

	if(month <= 2) {
		y = year - 1;
		m = month + 12;
	}
	else {
		y = year;
		m = month;
	}

	A = (long) (y / 100.);
	B = 2 - A + (long) (A / 4.);

	jd = (long) (365.25 * y) + (long) (30.6001 * (m + 1)) + day +
		1720994.5;

	jd += hour / 24. + minute / 1440. + second / 86400.;

	if(year > 1583) return(jd + B);
		else return(jd);

/*
	int jd1=367*year-7*(year+(month+9)/12)/4+275*month/9+day;
	double jd=jd1+1721013.5;
	jd+=(hour+minute/60.0+second/3600.0)/86400.0;
	
	return jd;
*/

}

/*
	JD2YMD  --  Convert	Julian	date  to  year,  month, day, 
	which are returned via integer pointers to integers (note that year is a long).
*/
void CDateTime::JD2YMD(double jd, int *yy, int *mm, int *dd)
{
	double z, f, a, alpha, b, c, d, e;

	jd += 0.5;
	z = floor(jd);
	f = jd - z;

	if (z < 2299161.0) {
		a = z;
	} else {
		alpha = floor((z - 1867216.25) / 36524.25);
		a = z + 1 + alpha - floor(alpha / 4);
	}

	b = a + 1524;
	c = floor((b - 122.1) / 365.25);
	d = floor(365.25 * c);
	e = floor((b - d) / 30.6001);

	*dd = (int) (b - d - floor(30.6001 * e) + f);
	*mm = (int) ((e < 14) ? (e - 1) : (e - 13));
	*yy = (long) ((*mm > 2) ? (c - 4716) : (c - 4715));
}

/*
	JHMS  --  Convert Julian time to hour, minutes, and seconds. 
*/
void CDateTime::JD2HMS(double jd, int *h, int *m, int *s)
{
    long ij;

    jd += 0.5;			      /* Astronomical to civil */
    ij = (long) (((jd - floor(jd)) * 86400.0) + 0.5);  // Round to nearest second
    *h = (int) (ij / 3600L);
    *m = (int) ((ij / 60L) % 60L);
    *s = (int) (ij % 60L);
}

/*
Returns true if y a leap year, and false otherwise, according
to the Gregorian calendar
*/
bool CDateTime::IsLeap(int y) 
{
	bool a = false;  
    if(y % 4 == 0) a = true;
    if(y % 100 == 0) a = false;      
    if(y % 400 == 0) a = true;
    return(a);
}

/* 
	returns day of week for a jd, 0 = Mon, 6 = Sun. 
*/
short CDateTime::Day_of_Week(double jd)
{
	double x;
	long i;
	short d;

	jd = jd+0.5;
	i =(long) jd; /* truncate */
	x = i/7.+0.01;
	d = (short)(7.*(x - (long) x));   /* truncate */
	return(d);
}

/*
	Adjust the input time to its normal range
*/
void CDateTime::AdjustTime(int *year, int *month, int *day, int *hour, int *minute, int *second)
{
	int mday[]={31,28,31,30,31,30,31,31,30,31,30,31};

	if((*second)>=59.999999999)
	{
		int quotient=((int)(*second))/60;
		(*second)-=60*quotient;
		(*minute)+=quotient;
	}
	if((*minute)>=59.999999999)
	{
		int quotient=((int)(*minute))/60;
		(*minute)-=60*quotient;
		(*hour)+=quotient;
	}
	if((*hour)>=24)
	{
		int quotient=((int)(*hour))/24;
		(*hour)-=24*quotient;
		(*day)+=quotient;
	}
	if(IsLeap(*year)) mday[1]=29;

	for(int i=11;i>=0;i--)
	{
		if((*month)==i+1)
		{
			if((*day)>mday[i]) 
			{
				int quotient=(*day)/mday[i];
				(*day)-=mday[i]*quotient;
				(*month)+=quotient;
				break;
			}
		}
	}

	if((*month)>12)
	{
		int quotient=(*month)/12;
		(*month)-=12*quotient;
		(*year)+=quotient;
		if((*month)==0)
		{
			(*month)=1;
		}
	}
}


//Special Operations begin////////////////////////////////////////////////////////////////
/*
	设置历元时间
*/
void CDateTime::SetEpochTime(short year, short month, short day, short hour, short minute, short second)
{
	m_epochjd=JD(year,month,day,hour,minute,second);
}

//设置历元时间
void CDateTime::SetEpochTime(double jd)
{
	m_epochjd=jd;
}

/*
	设置起始时间
*/
void CDateTime::SetStartTime(short year, short month, short day, short hour, short minute, short second)
{
	m_startjd=JD(year,month,day,hour,minute,second);
	m_currentjd=m_startjd;
}

//设置起始时间
void CDateTime::SetStartTime(double jd)
{
	m_startjd=jd;
	m_currentjd=m_startjd;
}

/*
	获取历元时间
*/
double CDateTime::GetEpochTime(BOOL IsSimMode)
{
	if(IsSimMode)
		return m_epochjd;
	else
		return CurrentTime(IsSimMode);
}

/*
	获取起始时间
*/
double CDateTime::GetStartTime(BOOL IsSimMode)
{
	if(IsSimMode)
		return m_startjd;
	else
		return CurrentTime(IsSimMode);
}

/*
	设置结束时间
*/
void CDateTime::SetEndTime(short year, short month, short day, short hour, short minute, short second)
{
	m_endjd=JD(year,month,day,hour,minute,second);
}

//设置结束时间
void CDateTime::SetEndTime(double jd)
{
	m_endjd=jd;
}

/*
	获取结束时间
*/
double CDateTime::GetEndTime(BOOL IsSimMode)
{
	if(IsSimMode)
		return m_endjd;
	else
		return CurrentTime(IsSimMode)+1;
}

/*
	获取当前时间
*/
double CDateTime::CurrentTime(BOOL IsSimMode)
{
	if(!IsSimMode)
	{
		SYSTEMTIME ut;
		GetSystemTime(&ut);
		m_year=ut.wYear ;
		m_month=ut.wMonth ;
		m_day=ut.wDay ;
		m_hour=ut.wHour ;
		m_minute=ut.wMinute ;
		m_second=ut.wSecond ;

		return JD(m_year,m_month,m_day,m_hour,m_minute,m_second);
	}
		return m_currentjd;
}

void CDateTime::SetStep(double step)
{
	m_stepjd=step;	
}

double CDateTime::GetStep()
{
	return m_stepjd;
}

void CDateTime::StepIt(BOOL IsSimMode)
{

	if(!IsSimMode)
		return;

	m_currentjd+=m_stepjd;

}

void CDateTime::Slowdown()
{
	m_stepjd/=2.;
}

void CDateTime::Speedup()
{
	m_stepjd*=2.;
}

void CDateTime::ResetTime()
{
	m_currentjd=m_startjd;
}
