#ifndef __FIR_H__
#define __FIR_H__
#include "DFD.h"
#include <vector>

enum enWindType
{
	WINDOW_RECT		= 0,
	WINDOW_BARTLETT = 1,
	WINDOW_HANNING	= 2,
	WINDOW_HAMMING	= 3,
	WINDOW_BLACKMAN = 4,
	WINDOW_GAUSS	= 5,
	WINDOW_KAISER	= 6,
};

const int		MAXTERM = 20;					//用于kaiser滤波器

/*************************窗函数声明**************************************/

template<typename Type> std::vector<Type> window(const enWindType&, int, Type);
template<typename Type> std::vector<Type> window(const enWindType&, int, Type, Type);

template<typename Type> std::vector<Type> kaiser(int, Type, Type);
template<typename Type> std::vector<Type> gauss(int, Type, Type);

template<typename Type> Type I0(Type alpha);



class CFIR : public CDFD
{
public:
	CFIR(const enFiltType &select, const enWindType &win);
	CFIR(const enFiltType &select, const enWindType &win, double a);
	~CFIR();

public:
	void    design();						//设计、生成滤波器(系数)
	void    dispInfo() const;				//输出显示设计的滤波器结果
	std::vector<double> getCoefs();			//获取滤波器的系数

private:
	void    orderEst();					//估算滤波器的长度(阶数)  注：严格的说FIR滤波器的阶数=其系数长度-1，此处认为相同不影响
	void    idealCoef();				//计算理想FIR滤波器系数
	void    calcCoef();					//计算实际FIR滤波器系数(理想的加窗)
	double  frqeResp(double freq);		//计算频率响应；计算离散时间傅里叶变换
	void    calcGain();					//计算增益
	bool    isSatisfy();				//判断设计的滤波器是否符合指标

	enWindType  m_windType;				// 窗口类型
	std::vector<double>  m_wind;		// window function
	std::vector<double>	m_coefs;		// coefficients
	std::vector<double>	m_edgeGain;		// 边界频率的增益，用来判断设计的滤波器是否符合指标
	double  m_alpha;					// window parameter

};	// class CFIR



#endif	//__FIR_H__


