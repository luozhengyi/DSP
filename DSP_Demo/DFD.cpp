#include "stdafx.h"
#include "DFD.h"
#include <iostream>

CDFD::CDFD(const enFiltType &select) : m_filtType(select)
{
}

CDFD::~CDFD()
{
}

double CDFD::getValue(double x, double min, double max)
{

	if ((x > min) && (x < max))
		return x;
	else
	{
		std::cout << "The parameter is out of range!" << std::endl;
		return 0;
	}
}

void CDFD::setParams(double fs, double f1, double a1, double f2, double a2)
{
	fsamp = getValue(fs, FREQ_MIN, FREQ_MAX);
	double maxFreq = fsamp / 2;		//根据奈奎斯特定理，信号的最大频率只能为采样频率的1/2

	if (m_filtType == FILTER_LOWPASS)
	{
		wpass1 = getValue(f1, FREQ_MIN, maxFreq);
		apass1 = getValue(a1, GAIN_TRAN, GAIN_PASS);
		wstop1 = getValue(f2, wpass1, maxFreq);
		astop1 = getValue(a2, GAIN_STOP, GAIN_TRAN);
	}
	else if (m_filtType == FILTER_HIGHPASS)
	{
		wstop1 = getValue(f1, FREQ_MIN, maxFreq);
		astop1 = getValue(a1, GAIN_STOP, GAIN_TRAN);
		wpass1 = getValue(f2, wstop1, maxFreq);
		apass1 = getValue(a2, GAIN_TRAN, GAIN_PASS);
	}
	else
		std::cout << "Parameters setting has failed!" << std::endl;

	// 将模拟的边界频率都转化为角频率(数字频率)
	wpass1 *= TWOPI;
	wstop1 *= TWOPI;

}

void CDFD::setParams(double fs, double f1, double a1, double f2,
	double f3, double a2, double f4, double a3)
{
	fsamp = getValue(fs, FREQ_MIN, FREQ_MAX);
	double maxFreq = fsamp / 2;

	if (m_filtType == FILTER_BANDPASS)
	{
		wstop1 = getValue(f1, FREQ_MIN, maxFreq);
		astop1 = getValue(a1, GAIN_STOP, GAIN_TRAN);
		wpass1 = getValue(f2, wstop1, maxFreq);
		wpass2 = getValue(f3, wpass1, maxFreq);
		apass1 = getValue(a2, GAIN_TRAN, GAIN_PASS);
		wstop2 = getValue(f4, wpass2, maxFreq);
		astop2 = getValue(a3, GAIN_STOP, GAIN_TRAN);
		apass2 = apass1;
		if (astop1 < astop2)
			astop2 = astop1;
		else
			astop1 = astop2;
	}
	else if (m_filtType == FILTER_BANDSTOP)
	{
		wpass1 = getValue(f1, FREQ_MIN, maxFreq);
		apass1 = getValue(a1, GAIN_TRAN, GAIN_PASS);
		wstop1 = getValue(f2, wpass1, maxFreq);
		wstop2 = getValue(f3, wstop1, maxFreq);
		astop1 = getValue(a2, GAIN_STOP, GAIN_TRAN);
		wpass2 = getValue(f4, wstop2, maxFreq);
		apass2 = getValue(a3, GAIN_TRAN, GAIN_PASS);
		astop2 = astop1;
		if (apass1 < apass2)
			apass1 = apass2;
		else
			apass2 = apass1;
	}
	else
		std::cout << "Parameters setting has failed!" << std::endl;

	// 将模拟的边界频率都转化为角频率(数字频率)
	wpass1 *= TWOPI;
	wpass2 *= TWOPI;
	wstop1 *= TWOPI;
	wstop2 *= TWOPI;
}
