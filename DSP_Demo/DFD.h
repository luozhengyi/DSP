#ifndef __DFD_H__
#define __DFD_H__
#include <string>

const double    FREQ_MIN = 1.0E-3;      // 最小频率
const double    FREQ_MAX = 1.0E12;      // 最大频率
const double    GAIN_PASS = -1.0E-2;    // 最大通带增益
const double    GAIN_TRAN = -3.0103;    // min pb, max sb gain
const double    GAIN_STOP = -1.0E02;    // 最小阻带增益

const double	EPS = 2.220446049250313e-016;	//浮点数最小的间隔，也即是浮点数的精度

const double	PI = 3.141592653589793;			//圆周率pi
const double	TWOPI = 6.283185307179586;


enum enFiltType
{
	FILTER_LOWPASS	= 0,
	FILTER_HIGHPASS = 1,
	FILTER_BANDPASS = 2,
	FILTER_BANDSTOP = 3,
};


class CDFD
{
public:
	CDFD(const enFiltType &select);
	~CDFD();
public:
	//设置滤波器参数：适用于低通和高通
	void setParams(double fs, double f1, double a1,double f2, double a2);
	//设置滤波器参数：适用于带通和带阻
	void setParams(double fs, double f1, double a1, double f2,double f3, double a2, double f4, double a3);

	//设计滤波器
	virtual void design() = 0;
	virtual void dispInfo() const = 0;

protected:

	//单位HZ
	double  fsamp;              // 采样频率
	double  wpass1, wpass2;     // 通带边界频率
	double  wstop1, wstop2;     // 阻带边界频率

	//单位db
	double  apass1, apass2;     // 通带增益
	double  astop1, astop2;     // 阻带增益

	enFiltType  m_filtType;     // 滤波器的类型
	int     order;              // 滤波器的长度

	//获取传入的参数
	inline double getValue(double x, double min, double max);

};

#endif	//__DFD_H__


