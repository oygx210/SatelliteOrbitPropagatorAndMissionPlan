// DateTime.h: interface for the CDateTime class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DATETIME_H_INCLUDED_)
#define AFX_DATETIME_H_INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "./Include/Constant.h"

#define BOOL bool

class CDateTime  
{
public:
	CDateTime();
	virtual ~CDateTime();

//Common Operations
public:
	double JD(long year, int month, int day, int hour, int minute, int second);
	void JD2YMD(double jd, int *yy, int *mm, int *dd);
	void JD2HMS(double jd, int *h, int *m, int *s);
	bool IsLeap(int y);
	short Day_of_Week(double jd);
	void AdjustTime(int *year, int *month, int *day, int *hour, int *minute, int *second);

//Special Operations
public:
	void SetEpochTime(short year, short month, short day, short hour, short minute, short second);
	void SetEpochTime(double jd);
	void SetStartTime(short year, short month, short day, short hour, short minute, short second);
	void SetStartTime(double jd);
	void SetEndTime(short year, short month, short day, short hour, short minute, short second);
	void SetEndTime(double jd);
	void SetStep(double step);
	void StepIt(BOOL IsSimMode=true);
	void Slowdown();
	void Speedup();
	double GetEpochTime(BOOL IsSimMode= true);
	double GetStartTime(BOOL IsSimMode= true);
	double CurrentTime(BOOL IsSimMode= true);
	double GetEndTime(BOOL IsSimMode= true);
	double GetStep();
	void ResetTime(void);


//Attributes
public:
	short m_year;
	short m_month;
	short m_day;
	short m_hour;
	short m_minute;
	short m_second;

	double m_epochjd;
	double m_startjd;
	double m_currentjd;
	double m_endjd;
	double m_stepjd;

};

#endif
