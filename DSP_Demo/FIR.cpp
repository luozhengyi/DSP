#include "stdafx.h"
#include "FIR.h"
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <math.h>
using namespace std;

using std::setw;
using std::ios;
using std::setiosflags;
using std::setprecision;


CFIR::CFIR(const enFiltType &select, const enWindType &win)
	: CDFD(select), m_windType(win)
{
	bool cond = (m_filtType >= FILTER_LOWPASS || m_filtType <= FILTER_BANDSTOP);
	assert(cond);

	cond = (m_windType >= WINDOW_RECT || m_windType <= WINDOW_GAUSS);
	assert(cond);
}

CFIR::CFIR(const enFiltType &select, const enWindType &win, double a)
	: CDFD(select), m_windType(win), m_alpha(a)
{
	bool cond = (m_filtType >= FILTER_LOWPASS || m_filtType <= FILTER_BANDSTOP);
	assert(cond);

	cond = (m_windType == WINDOW_KAISER || m_windType == WINDOW_GAUSS);
	assert(cond);
}

CFIR::~CFIR()
{
}



void CFIR::design()
{
	// Estimate the length of filter.
	orderEst();		//估算滤波器的长度，阶数
	calcCoef();		//计算滤波器的系数
	calcGain();		//计算增益
	while (!isSatisfy())
	{
		order += 4;
		calcCoef();
		calcGain();
	}
		
}

void CFIR::dispInfo() const
{
	int i, j, k;

	cout << std::endl;
	cout << "\t\t    Filter selectivity      :  " << m_filtType << endl;
	cout << "\t\t    Window type             :  " << m_windType << endl;
	cout << "\t\t    Sampling Frequency (Hz) :  " << fsamp << endl;

	//  gains and edge frequency
	if (m_filtType == FILTER_LOWPASS)
	{
		cout << "\t\t    Passband frequency (Hz) :  "
			<< wpass1 / TWOPI << endl;
		cout << "\t\t    Passband gain      (dB) :  "
			<< apass1 << endl;
		cout << "\t\t    Stopband frequency (Hz) :  "
			<< wstop1 / TWOPI << endl;
		cout << "\t\t    Stopband gain      (dB) :  "
			<< astop1 << endl;
	}
	else if (m_filtType == FILTER_HIGHPASS)
	{
		cout << "\t\t    Stopband frequency (Hz) :  "
			<< wstop1 / TWOPI << endl;
		cout << "\t\t    Stopband gain      (dB) :  "
			<< astop1 << endl;
		cout << "\t\t    Passband frequency (Hz) :  "
			<< wpass1 / TWOPI << endl;
		cout << "\t\t    Passband gain      (dB) :  "
			<< apass1 << endl;
	}
	else if (m_filtType == FILTER_BANDPASS)
	{
		cout << "\t\t    Lower stopband frequency (Hz) :  "
			<< wstop1 / TWOPI << endl;
		cout << "\t\t    Lower stopband gain      (dB) :  "
			<< astop1 << endl;
		cout << "\t\t    Lower passband frequency (Hz) :  "
			<< wpass1 / TWOPI << endl;
		cout << "\t\t    Upper passband frequency (Hz) :  "
			<< wpass2 / TWOPI << endl;
		cout << "\t\t    Passband gain            (dB) :  "
			<< apass1 << endl;
		cout << "\t\t    Upper stopband frequency (Hz) :  "
			<< wstop2 / TWOPI << endl;
		cout << "\t\t    Upper stopband gain      (dB) :  "
			<< astop2 << endl;
	}
	else
	{
		cout << "\t\t    Lower passband frequency (Hz) :  "
			<< wpass1 / TWOPI << endl;
		cout << "\t\t    Lower passband gain      (dB) :  "
			<< apass1 << endl;
		cout << "\t\t    Lower stopband frequency (Hz) :  "
			<< wstop1 / TWOPI << endl;
		cout << "\t\t    Upper stopband frequency (Hz) :  "
			<< wstop2 / TWOPI << endl;
		cout << "\t\t    Stopband gain            (dB) :  "
			<< astop1 << endl;
		cout << "\t\t    Upper passband frequency (Hz) :  "
			<< wpass2 / TWOPI << endl;
		cout << "\t\t    Upper passband gain      (dB) :  "
			<< apass2 << endl;
	}

	// display coefficients
	cout << endl << endl;
	cout << "\t\t\t\t  Filter Coefficients" << endl << endl;
	cout << "   N    [     N + 0            N + 1";
	cout << "            N + 2            N + 3     ]" << endl;
	cout << "  ===   ============================";
	cout << "========================================";
	for (i = 0; i<order / 4; ++i)
	{
		j = i * 4;
		cout << endl << setw(4) << j << "    ";

		cout << setiosflags(ios::scientific) << setprecision(8);
		for (k = 0; k<4; ++k)
			cout << setw(16) << m_coefs[j + k] << " ";
	}

	cout << endl << endl << endl
		<< "\t ==================== Edge Frequency Response";
	cout << " ====================" << endl;
	cout << setiosflags(ios::fixed);
	if (m_edgeGain.size() == 2)
	{
		cout << "\t     Mag(fp) = " << m_edgeGain[0] << "(dB)";
		cout << "       Mag(fs) = " << m_edgeGain[1] << "(dB)" << endl;
	}
	else
	{
		cout << "\t     Mag(fp1) = " << m_edgeGain[0] << "(dB)";
		cout << "       Mag(fp2) = " << m_edgeGain[1] << "(dB)" << endl;
		cout << "\t     Mag(fs1) = " << m_edgeGain[2] << "(dB)";
		cout << "       Mag(fs2) = " << m_edgeGain[3] << "(dB)" << endl;
	}
}

std::vector<double> CFIR::getCoefs()
{
	return m_coefs;
}

void CFIR::orderEst()
{

	double  deltaFreq,  // 过渡带宽
			lowerDF,    // 低过渡带宽(带通/带阻有两个过渡带)
			upperDF,    // 高过渡带宽
			errSB,      // stopband error
			errPB,      // passband error
			errMin,     // minimum of sb/pb errors
			errDB;      // minimum error in dB

	// Determine frequency delta
	if (m_filtType == FILTER_LOWPASS || m_filtType == FILTER_HIGHPASS)
		deltaFreq = abs(wstop1 - wpass1) / fsamp;
	else
	{
		lowerDF = abs(wstop1 - wpass1) / fsamp;
		upperDF = abs(wstop2 - wpass2) / fsamp;

		if (lowerDF > upperDF)
			deltaFreq = upperDF;
		else
			deltaFreq = lowerDF;
	}

	// Determine stopband and passband errors
	errPB = 1 - pow(10, 0.05*apass1);		//apass1=20*lg(1-errPB)；errPB:通带衰减；apass1:通带波纹幅度
	errSB = pow(10, 0.05*astop1);			//astop1=20*lg(errSB)：阻带衰减

	if (errSB < errPB)
		errMin = errSB;
	else
		errMin = errPB;
	errDB = -20 * log10(errMin);

	// Store filter length in pFilt and return beta.
	if (errDB > 21)
		order = int(ceil(1 + (errDB - 7.95) / (2.285*deltaFreq)));
	else
		order = int(ceil(1 + (5.794 / deltaFreq)));

	order += 4 - order % 4;
	//order = 2 * ( order / 2) + 1;	//滤波器阶数为奇数
}

void CFIR::idealCoef()
{
	int     i;
	double  t,
			tau,		//延迟
			Wc1,		//过渡带中心频率(带通或带阻有两个)
			Wc2;		//过渡带中心频率

	// Calculate tau as non-integer if order is even.
	tau = (order - 1) / 2.0;

	//将截止频率设为通带边界频率和阻带边界频率的1/2;同时将其转化为数字频率(归一化频率：0到1，或者0到pi)
	Wc1 = (wstop1 + wpass1) / (2 * fsamp);
	Wc2 = (wstop2 + wpass2) / (2 * fsamp);

	// Calc coefs based on selectivity of filter.
	if (m_filtType == FILTER_LOWPASS)
		for (i = 0; i<order; ++i)
		{
			if (i == tau)
				m_coefs[i] = Wc1 / PI;
			else
			{
				t = i - tau;
				m_coefs[i] = sin(Wc1*t) / (PI*t);
			}
		}
	else if (m_filtType == FILTER_HIGHPASS)
		for (i = 0; i<order; ++i)
		{
			if (i == tau)
				m_coefs[i] = (PI - Wc1) / PI;
			else
			{
				t = i - tau;
				m_coefs[i] = (sin(PI*t) - sin(Wc1*t)) / (PI*t);
			}
		}
	else if (m_filtType == FILTER_BANDPASS)
		for (i = 0; i<order; ++i)
		{
			if (i == tau)
				m_coefs[i] = (Wc2 - Wc1) / PI;
			else
			{
				t = i - tau;
				m_coefs[i] = (sin(Wc2*t) - sin(Wc1*t)) / (PI*t);
			}
		}
	else
		for (i = 0; i<order; ++i)
		{
			if (i == tau)
				m_coefs[i] = (PI + Wc1 - Wc2) / PI;
			else
			{
				t = i - tau;
				m_coefs[i] = (sin(PI*t) - sin(Wc2*t) + sin(Wc1*t)) / (PI*t);
			}
		}
}

void CFIR::calcCoef()
{
	m_coefs.resize(order);
	m_wind.resize(order);

	//  Calculate the ideal FIR coefficients.
	idealCoef();

	// Get window function.
	if (m_windType == WINDOW_KAISER || m_windType == WINDOW_GAUSS)
		m_wind = window(m_windType, order, m_alpha, 1.0);
	else
		m_wind = window(m_windType, order, 1.0);

	//  Multiply window and ideal coefficients.
	for (unsigned i = 0; i < m_coefs.size(); i++)
	{
		m_coefs[i] *= m_wind[i];
	}
	
}

double CFIR::frqeResp(double freq)
{
	double  mag = 1.0,
			rea = 0.0,
			img = 0.0,
			omega = TWOPI*freq / fsamp;		//角频率=2*pi*f/fs

	// 计算离散时间傅里叶变换 当（w=omega）时 X(ejw)的值
	for (int i = 0; i<order; ++i)
	{
		double  domega = i * omega;
		rea += m_coefs[i] * cos(domega);
		img += m_coefs[i] * sin(domega);
	}

	// Calculate final result and convert to degrees.
	mag = sqrt(rea*rea + img*img);
	if (mag < EPS)
		mag = EPS;
	mag = 20 * log10(mag);

	return mag;
}

void CFIR::calcGain()
{
	// Determine the edge freqency.
	if (m_filtType == FILTER_LOWPASS || m_filtType == FILTER_HIGHPASS)
	{
		double  f1 = wpass1 / TWOPI,
			f2 = wstop1 / TWOPI;

		m_edgeGain.resize(2);
		m_edgeGain[0] = frqeResp(f1);
		m_edgeGain[1] = frqeResp(f2);
	}
	else
	{
		double  f1 = wpass1 / TWOPI,
				f2 = wpass2 / TWOPI,
				f3 = wstop1 / TWOPI,
				f4 = wstop2 / TWOPI;

		m_edgeGain.resize(4);
		m_edgeGain[0] = frqeResp(f1);
		m_edgeGain[1] = frqeResp(f2);
		m_edgeGain[2] = frqeResp(f3);
		m_edgeGain[3] = frqeResp(f4);
	}
}

bool CFIR::isSatisfy()
{
	if (m_edgeGain.size() == 2)
	{
		if (m_edgeGain[0]<apass1 || m_edgeGain[1]>astop1)
			return false;
	}
	else
	{
		if (m_edgeGain[0]<apass1 || m_edgeGain[1]<apass2 ||
			m_edgeGain[2]>astop1 || m_edgeGain[3]>astop2)
			return false;
	}

	return true;
}





/*****************************窗函数设计*********************************/

template <typename Type>
std::vector<Type> window(const enWindType &wnName, int N, Type amp)
{
	std::vector<Type> win(N);

	if (wnName == WINDOW_RECT)
		for (int i = 0; i < (N + 1) / 2; ++i)
		{
			win[i] = amp;
			win[N - 1 - i] = win[i];
		}
	else if (wnName == WINDOW_BARTLETT)
		for (int i = 0; i < (N + 1) / 2; ++i)
		{
			win[i] = amp * 2 * i / (N - 1);
			win[N - 1 - i] = win[i];
		}
	else if (wnName == WINDOW_HANNING)
		for (int i = 0; i < (N + 1) / 2; ++i)
		{
			win[i] = amp * Type(0.5 - 0.5*cos(TWOPI*i / (N - 1)));
			win[N - 1 - i] = win[i];
		}
	else if (wnName == WINDOW_HAMMING)
		for (int i = 0; i < (N + 1) / 2; ++i)
		{
			win[i] = amp * Type(0.54 - 0.46*cos(TWOPI*i / (N - 1.0)));
			win[N - 1 - i] = win[i];
		}
	else if (wnName == WINDOW_BLACKMAN)
		for (int i = 0; i < (N + 1) / 2; ++i)
		{
			win[i] = amp * Type(0.42 - 0.50*cos(TWOPI*i / (N - 1.0)) + 0.08*cos(2 * TWOPI*i / (N - 1.0)));
			win[N - 1 - i] = win[i];
		}
	else
	{
		cerr << "No such type window!" << endl;
		return std::vector<Type>(0);
	}

	return win;
}

template <typename Type>
std::vector<Type> window(const enWindType &wnName, int N, Type alpha, Type amp)
{
	if (wnName == WINDOW_KAISER)
		return kaiser(N, alpha, amp);
	else if (wnName == WINDOW_GAUSS)
		return gauss(N, alpha, amp);
	else
	{
		cerr << "No such type window!" << endl;
		return std::vector<Type>(0);
	}
}


/**
* Calculates kaiser window coefficients.
*/
template <typename Type>
std::vector<Type> kaiser(int N, Type alpha, Type amp)
{
	std::vector<Type> win(N);

	for (int i = 0; i < (N + 1) / 2; ++i)
	{
		Type beta = 2 * alpha * Type(sqrt(i*(N - i - 1.0)) / (N - 1.0));
		win[i] = amp * I0(beta) / I0(alpha);
		win[N - 1 - i] = win[i];
	}

	return win;
}


/**
* Calculates gauss window coefficients. "Alpha: is a optional parameter,
* the default value is 2.5.
*/
template <typename Type>
std::vector<Type> gauss(int N, Type alpha, Type amp)
{
	std::vector<Type> win(N);
	Type center = (N - 1) / Type(2);

	for (int i = 0; i < (N + 1) / 2; ++i)
	{
		Type tmp = alpha*(i - center) / center;
		win[i] = amp * Type(exp(-0.5*tmp*tmp));
		win[N - 1 - i] = win[i];
	}

	return win;
}


/**
* The zeroth N modified Bessel function of the first kind.
*/
template <typename Type>
Type I0(Type alpha)		//是i零，不是io
{
	double  J = 1.0,
		K = alpha / 2.0,
		iOld = 1.0,
		iNew;
	bool    converge = false;

	// Use series expansion definition of Bessel.
	for (int i = 1; i < MAXTERM; ++i)
	{
		J *= K / i;
		iNew = iOld + J*J;

		if ((iNew - iOld) < EPS)
		{
			converge = true;
			break;
		}
		iOld = iNew;
	}

	if (!converge)
		return Type(0);

	return Type(iNew);
}