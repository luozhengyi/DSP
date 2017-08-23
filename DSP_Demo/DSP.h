#ifndef __DSP_H__
#define __DSP_H__
#include <vector>
#include "FIR.h"
//fftw引用
#include "fftw\\fftw3.h"
#pragma comment(lib, "fftw\\libfftw3-3.lib") // double版本

class CDSP
{
public:
	CDSP();
	~CDSP();

public:
	/**
	* @fn	bool QAM_Gen();
	* @brief	生成16QAM信号.
	* @param	pIData:生成的I路数据	pQData:生成的Q路数据	fs:采样率	fd:符号率	numSymbols:符号数
	* @return	true if it succeeds,false if it fails
	*/
	bool QAM_Gen(double fs, double fd, double numSymbols);

	/**
	* @fn	bool QAM_Demod();
	* @brief	解调16QAM信号.
	* @param	pData:I、Q两路合成信号	iDataLen:数据长度	fs:采样率	fc:载波频率		fd:符号率	R:滚降系数
	* @return	
	*/
	bool QAM_Demod(double *pData, int iDataLen, double fs, double fc, double fd, double R=0.5);
public:
	//xn的长度要大于等于hn的长度！
	std::vector<double> Conv(double* xn, int xn_len, double* hn, int hn_len);		//利用卷积定义求线性卷积
	std::vector<double> FastConv(double* xn, int xn_len, double* hn, int hn_len);	//利用fft来求循环(线性)卷积
public:
	void FilterFilt(double *pData, int iDataLen, std::vector<double>& filter);		//执行滤波,类似于matlab中filter,取前n个数，会造成相偏
	void FiltFilt(double *pData, int iDataLen, std::vector<double>& filter);		//等价于matlab中的filtfilt()零相偏滤波器
public:
	void SaveToFile(double *pData, int iDataLen);					//将数据保存到文件
	void SaveToFile(double *pIData, double *pQData, int iDataLen);	//将数据保存到文件

public:
	std::vector<double> RCosine(double fs, double fd, double R);	//升余弦滤波器
public:
	void FilterIR(std::vector<double> coefs);	//滤波器单位冲击响应
public:
	std::vector<double> CDSP::GetDataFromFile();	//获取数据
	void DemDataFromFile();
	


private:
	double *m_pIData;
	double *m_pQData;
};
#endif	//__DSP_H__


