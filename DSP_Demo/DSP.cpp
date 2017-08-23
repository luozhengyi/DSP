#include "stdafx.h"
#include "DSP.h"
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <algorithm>

CDSP::CDSP()
{
	m_pIData = NULL;
	m_pQData = NULL;
}

CDSP::~CDSP()
{
	if (m_pIData != NULL)
		delete[] m_pIData;
	if (m_pQData != NULL)
		delete[] m_pQData;
}

bool CDSP::QAM_Gen(double fs, double fd, double numSymbols)
{
	int oversampling = fs / fd;					//过采样率
	int iDataLen = numSymbols * oversampling;	//I、Q路数据长度
	int Ref = 1;								//参考电压
	double R = 0.5;								//滚降系数

	//分配内存
	double *pIData = new double[iDataLen]();	//将分配的内存初始化为0
	double *pQData = new double[iDataLen]();
	if (pIData == NULL || pQData == NULL)
		return false;
	memset(pIData, 0.0, iDataLen*sizeof(double));
	memset(pQData, 0.0, iDataLen*sizeof(double));

	//产生随机数（完成插值）
	srand((unsigned)time(NULL));	//初始化随机数种子
	for (int i = 0, index = (i * oversampling); index < iDataLen; i++, index = (i * oversampling))
	{
		pIData[index] = rand() % 4 - 1.5;	//产生-1.5 -0.5 0.5 1.5的随机数
		pQData[index] = rand() % 4 - 1.5;	//产生-1.5 -0.5 0.5 1.5的随机数
	}
	

	//低通滤波：截止频率(fd =1MHZ)（截止频率太大太小都不一定是好事,经过测试）
	double  fpass = (fd - 0.2e6) - 0.05e6,
			apass = -3,
			fstop = (fd - 0.2e6) + 0.05e6,
			astop = -20;		//-10db要比-20db滤波效果好，但是对信号却并不一定是好事
	CFIR fir(FILTER_LOWPASS, WINDOW_BLACKMAN);
	fir.setParams( fs, fpass, apass, fstop, astop );
	fir.design();
	//fir.dispInfo();
	std::vector<double>	coefs = fir.getCoefs();
	//FilterIR(coefs);		//单位冲击响应
	FilterFilt(pIData, iDataLen, coefs);
	FilterFilt(pQData, iDataLen, coefs);

	//保存数据到文件
	SaveToFile(pIData, pQData, iDataLen);

	//基带脉冲成型
	coefs = RCosine(fs, fd, R);
	FilterIR(coefs);
	FilterFilt(pIData, iDataLen, coefs);
	FilterFilt(pQData, iDataLen, coefs);
	//保存数据到文件
	SaveToFile(pIData, pQData, iDataLen);

	//将基带信号调制到高频上
	double fc = 20e6;		//功率一旦太大会有谐波产生！
	for (int i = 0; i < iDataLen; i++)
	{
		/*pIData[i] *= cos(TWOPI*fc / fs*i);
		pQData[i] *= -sin(TWOPI*fc / fs*i);*/

		//matlab下标从1开始,这其实就是时延，就会造成相偏
		pIData[i] *= cos(TWOPI*fc / fs*(i+1));
		pQData[i] *= -sin(TWOPI*fc / fs*(i + 1));

		pIData[i] += pQData[i];	//将I、Q路结果叠加到一路存放到I路的内存上。
	}
	

	//将数据写入文件
	SaveToFile(pIData, iDataLen);

	//解调
	QAM_Demod(pIData, iDataLen, fs, fc, fd, R);

	
	
	//滤波器测试
	//FilterIR(coefs);


	

	//释放资源
	delete[] pIData;
	delete[] pQData;
	return true;
}

bool CDSP::QAM_Demod(double *pData,int iDataLen, double fs, double fc, double fd, double R)
{
	int oversampling = fs / fd;					//过采样率
	int numSymbols = iDataLen / oversampling;	//符号数

	//开辟内存
	double *pIData = new double[iDataLen]();
	double *pQData = new double[iDataLen]();
	double *pISymbols = new double[numSymbols]();
	double *pQSymbols = new double[numSymbols]();
	if (pIData == NULL || pQData == NULL || pISymbols == NULL || pQSymbols == NULL)
		return false;

	//带通滤波，滤出载波附近的信号;带通滤波不能滤得太狠，否则会发生相偏,等各种问题(FilterFilt()会造成相偏)
	double  fstop1 = fc - fd * (1.0 + R) - 0.1e6,
			astop1 = -20,
			fpass1 = fc - fd * (1.0 + R) + 0.1e6,
			fpass2 = fc + fd * (1.0 + R) - 0.1e6,
			apass1 = -3,
			fstop2 = fc + fd * (1.0 + R) + 0.1e6,
			astop2 = -20;
	CFIR fir(FILTER_BANDPASS, WINDOW_BLACKMAN);
	fir.setParams(fs, fstop1, astop1, fpass1, fpass2, apass1, fstop2, astop2);
	fir.design();
	//fir.dispInfo();
	std::vector<double>	coefs = fir.getCoefs();
	//FilterFilt(pData, iDataLen, coefs);
	FiltFilt(pData, iDataLen, coefs);	//零相偏滤波
	//std::vector<double>	coefs;


	//相干解调，分离出I、Q两路信号
	for (int i = 0; i < iDataLen; i++)
	{
		/*pIData[i] = 2.0 * pData[i] * cos(TWOPI * fc / fs * i);
		pQData[i] = 2.0 * pData[i] * -sin(TWOPI * fc / fs * i);*/

		//matlab中下标是从1开始的,与matlab结果对比
		pIData[i] = 2.0 * pData[i] * cos(TWOPI * fc / fs * (i + 1));
		pQData[i] = 2.0 * pData[i] * -sin(TWOPI * fc / fs * (i + 1));
	}

	//保存数据到文件
	SaveToFile(pIData, pQData, iDataLen);

	//匹配滤波（也是低通滤波）
	coefs=RCosine(fs, fd, R);
	FilterFilt(pIData, iDataLen, coefs);
	FilterFilt(pQData, iDataLen, coefs);

	//保存数据到文件
	SaveToFile(pIData, pQData, iDataLen);

	

	//寻找最佳采样点,一个符号周期内采样了oversampling个点，以每个点作为抽样的起点，计算功率。
	int start = 0;		//采样点起始下标
	double *pPower = new double[oversampling]();	//保存功率
	if (pPower != NULL)
		for (int i = 0; i < oversampling; i++)
		{
			int index = 0;
			for (int j = 0; j < numSymbols/2; j++)		//可以不用取全部的符号
			{
				index = i + j*oversampling;
				pPower[i] += (pow(pIData[index],2) + pow(pQData[index],2));		//用绝对值代替平方是一样的
			}
				
		}
	start = std::max_element(pPower, pPower + oversampling ) - pPower;	//寻找功率最大值的索引,oversampling不需-1

	//根据最佳采样点进行抽样
	for (int i = 0, index = 0; i < numSymbols && index < iDataLen; ++i, index = start + i*oversampling)
	{
		pISymbols[i] = pIData[index];
		pQSymbols[i] = pQData[index];
	}


	//计算EVM


	//保存数据到文件
	SaveToFile(pISymbols, pQSymbols, numSymbols);

	//释放资源
	delete[] pIData;
	delete[] pQData;
	delete[] pISymbols;
	delete[] pQSymbols;
	delete[] pPower;
	return true;

}

std::vector<double> CDSP::FastConv(double* xn, int xn_len, double* hn, int hn_len)
{
	//xn_len >= hn_len
	if (xn_len < hn_len)
		return std::vector<double> (0);

	int iLen = xn_len + hn_len - 1;

	//重新分配内存
	//if (iLen % 2 != 0)	//fft变换需要是2的倍数(好像不需要哦)
	//iLen++;
	double *yn = (double *)fftw_malloc(sizeof(double) * iLen);
	std::vector<double> ve_yn(iLen, 0.0);
	double *xn_new = (double *)fftw_malloc(sizeof(double) * iLen);
	double *hn_new = (double *)fftw_malloc(sizeof(double) * iLen);
	fftw_complex *XK = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * iLen);
	fftw_complex *HK = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * iLen);
	fftw_complex *YK = NULL;

	//内存拷贝
	memmove(xn_new, xn, xn_len*sizeof(double));
	memmove(hn_new, hn, hn_len*sizeof(double));

	//进行fft变换ll
	fftw_plan p1, p2;
	p1 = fftw_plan_dft_r2c_1d(iLen, xn_new, XK, FFTW_ESTIMATE);	//实数到复数，由Hermite对称性只取fs的一半
	p2 = fftw_plan_dft_r2c_1d(iLen, hn_new, HK, FFTW_ESTIMATE);
	fftw_execute(p1);
	fftw_execute(p2);
	//根据实序列FFT的共轭对称性求FFT结果的后一半值
	for (int i = iLen / 2 + 1; i < iLen; i++)
	{
		XK[i][0] = XK[iLen - i][0];
		XK[i][1] = -XK[iLen - i][1];

		HK[i][0] = HK[iLen - i][0];
		HK[i][1] = -HK[iLen - i][1];
	}

	//进行ifft变换
	fftw_plan p3;
	YK = XK;
	double Re = 0.0;
	double Im = 0.0;
	for (int i = 0; i < iLen; i++)
	{
		Re = 0.0;
		Im = 0.0;
		Re = XK[i][0] * HK[i][0] - XK[i][1] * HK[i][1];
		Im = XK[i][0] * HK[i][1] + XK[i][1] * HK[i][0];
		YK[i][0] = Re / iLen;
		YK[i][1] = Im / iLen;
	}
	p3 = fftw_plan_dft_c2r_1d(iLen, YK, yn, FFTW_ESTIMATE);
	fftw_execute(p3);
	for (int i = 0; i < iLen; i++)	//将数组赋值给 vector
	{
		ve_yn[i] = yn[i];
	}

	// 释放资源
	fftw_destroy_plan(p1);
	fftw_destroy_plan(p2);
	fftw_destroy_plan(p3);
	fftw_free(HK);
	fftw_free(XK);
	fftw_free(xn_new);
	fftw_free(hn_new);
	fftw_free(yn);

	return ve_yn;
}

std::vector<double> CDSP::Conv(double* xn, int xn_len, double* hn, int hn_len)
{
	//xn_len >= hn_len
	if (xn_len < hn_len)
		return std::vector<double>(0);


	int iLen = xn_len + hn_len - 1;
	std::vector<double> yn(iLen,0.0);

	for (int i = 0; i < iLen; ++i)
	{
		yn[i] = 0;
		if (i <= hn_len-1)
			for (int j = 0; j <= i; ++j)
				yn[i] += hn[j] * xn[i - j];
		else if (i <= xn_len-1)
			for (int j = 0; j <= hn_len-1; ++j)
				yn[i] += hn[j] * xn[i - j];
		else
			for (int j = i - xn_len+1 ; j <= hn_len-1; ++j)
				yn[i] += hn[j] * xn[i - j];
	}
	return yn;
}

void CDSP::FilterFilt(double *pData, int iDataLen, std::vector<double>& filter)
{
	std::vector<double> result;

	//vector to Array
	double *pfilter = new double[filter.size()]();
	for (unsigned i = 0; i < filter.size(); i++)
		pfilter[i] = filter[i];
	
	//滤波(卷积）
	result = FastConv(pData, iDataLen, pfilter, filter.size());

	//滤波结果赋值给原始数据内存
	for (int i = 0; i < iDataLen; i++)
		pData[i] = result[i];
	//
	/*int index = 0;
	index = (filter.size()-1) / 2;
	for (int i = index; i < index+iDataLen; i++)
	{
		pData[i - index] = result[i];
	}*/
		

	//释放资源
	delete[] pfilter;
}

void CDSP::FiltFilt(double *pData, int iDataLen, std::vector<double>& filter)
{
	//判断数据是否有效
	int iFilterLen = filter.size();
	int iNewDataLen = 2*iFilterLen + iDataLen;	//不乘2似乎也可以
	if (iFilterLen <= 0 || pData  == NULL)
		return;

	//改变原始数据长度，分配内存
	double *pNewData = new double[iNewDataLen]();
	memcpy(pNewData, pData, iDataLen*sizeof(double));

	//第一次滤波(卷积）
	FilterFilt(pNewData,iNewDataLen,filter);
	//第一反转
	for (int i = 0; i < (iNewDataLen-1)/2; i++)
	{
		double temp = pNewData[i];
		pNewData[i] = pNewData[iNewDataLen - 1 - i];
		pNewData[iNewDataLen - 1 - i] = temp;
	}
	//第二次滤波
	FilterFilt(pNewData, iNewDataLen, filter);
	//第二反转
	for (int i = 0; i < (iNewDataLen - 1) / 2; i++)
	{
		double temp = pNewData[i];
		pNewData[i] = pNewData[iNewDataLen - 1 - i];
		pNewData[iNewDataLen - 1 - i] = temp;
	}


	//滤波结果赋值给原始数据内存
	memcpy(pData, pNewData, iDataLen*sizeof(double));

	//释放资源
	delete[] pNewData;
}

void CDSP::SaveToFile(double *pData, int iDataLen)
{
	if (pData == NULL)	//判断数据指针是否有效
		return;
	std::fstream file1;
	file1.open("E:\\iqtools_2015_02_07\\iqtools\\vc_DATA.bin", std::ios::out | std::ios::binary);
	if (!file1)
		return;

	double minData = *(std::min_element(pData, pData + iDataLen - 1));
	double maxData = *(std::max_element(pData, pData + iDataLen - 1));
	double scale = abs(maxData) > abs(minData) ? abs(maxData) : abs(minData);

	short data = 0;
	char ch_data[2];
	for (int i = 0; i < iDataLen; i++)
	{
		data = 0;
		memset(ch_data, 0, 2);
		//归一化会导致在matlab中画出来的幅度会有一个倍数关系
		//可以除以10降低20db
		data = (short)(pData[i] / scale * ((1 << 15) - 1));	//pData[i] / scale一定得是-1到1之间的范围，不然就溢出了。
		memcpy(ch_data, &data, 2);
		file1.write(ch_data, 2);
		
	}

	file1.close();
}

void CDSP::SaveToFile(double *pIData, double *pQData, int iDataLen)
{
	if (pIData == NULL || pQData == NULL)	//判断数据指针是否有效
		return;
	std::fstream file1;
	file1.open("C:\\Users\\machenike\\Desktop\\16QAM调制与解调\\Lab_M_file\\vc_IQDATA.bin", std::ios::out | std::ios::binary);
	if (!file1)
		return;

	double minIData = *(std::min_element(pIData, pIData + iDataLen - 1));
	double maxIData = *(std::max_element(pIData, pIData + iDataLen - 1));
	double i_scale = abs(maxIData) > abs(minIData) ? abs(maxIData) : abs(minIData);

	double minQData = *(std::min_element(pQData, pQData + iDataLen - 1));
	double maxQData = *(std::max_element(pQData, pQData + iDataLen - 1));
	double q_scale = abs(maxQData) > abs(minQData) ? abs(maxQData) : abs(minQData);

	short data = 0;
	char ch_data[2];
	for (int i = 0; i < iDataLen; i++)
	{
		data = 0;
		memset(ch_data, 0, 2);
		//归一化会导致在matlab中画出来的幅度会有一个倍数关系
		//可以除以10降低20db
		data = (short)(pIData[i] / i_scale * ((1 << 15) - 1));	//pData[i] / scale一定得是-1到1之间的范围，不然就溢出了。
		memcpy(ch_data, &data, 2);
		file1.write(ch_data, 2);

		data = 0;
		memset(ch_data, 0, 2);
		//归一化会导致在matlab中画出来的幅度会有一个倍数关系
		//可以除以10降低20db
		data = (short)(pQData[i] / q_scale * ((1 << 15) - 1));	//pData[i] / scale一定得是-1到1之间的范围，不然就溢出了。
		memcpy(ch_data, &data, 2);
		file1.write(ch_data, 2);

	}

	file1.close();
}

void CDSP::FilterIR(std::vector<double>	coefs)
{
	//冲击信号
	const int iDataLen = 40000;
	double dbData[iDataLen] = { 0.0 };
	dbData[0] = 1.0;

	FilterFilt(dbData, iDataLen, coefs);


	//保存至文件
	SaveToFile(dbData, iDataLen);
}

std::vector<double> CDSP::RCosine(double fs, double fd, double R)
{
	if (R <= 0 || R>1)		//验证滚降系数
		return std::vector<double>();

	int order = (int)6.0 * fs / fd +1 ;	//滤波器阶数（公式网上抄的）
	double tau = (order - 1) / 2.0;		//滤波器延迟
	double Wc1 = PI*fd / fs;			//截止频率
	double t=0.0;
	std::vector<double> coefs(order);
	for (int i = 0; i < order; ++i)
	{
		t = i - tau;
		if (t == 0)
			coefs[i] = (1.0 - R) + 4.0*R / PI;
		else if (t == 1.0 / (4.0 * fd / fs * R) || t == -1.0 / (4.0 * fd / fs * R))		//不要引入PI，精度有问题导致判断出错
			coefs[i] = R / sqrt(2.0) *((1 + 2 / PI)*sin(PI / (4.0 * R)) + (1 - 2 / PI)*cos(PI / (4.0 * R)));
		else
			coefs[i] = (4*R/PI*cos((1 + R)*Wc1*t) + sin((1 - R)*Wc1*t) / (Wc1*t)) / (1 - pow(4*R*t*Wc1/PI, 2));
	}

	return coefs;
}

std::vector<double> CDSP::GetDataFromFile()	//获取绘图数据
{
	std::vector<double> veData;
	char szDataPath[100] = "C:\\Users\\machenike\\Desktop\\cppData_transfer1.bin";
	int iDataLen = 0;
	do
	{
		if (strlen(szDataPath) == 0)
			break;
		std::fstream file1;
		file1.open(szDataPath, std::ios::in | std::ios::binary);
		if (!file1)
			break;
		//读取文件，写如文件，按行处理
		int num = 0;
		char ch[2];
		short sData = 0;
		while (!file1.eof())
		{
			memset(ch, 0, sizeof(ch));
			sData = 0;
			file1.read(ch, sizeof(ch));	//读取两个字节
			if (file1.good())			//防止到文件结束处多读取一次
			{
				memcpy(&sData, ch, 2);
				veData.push_back(sData);
				num++;
				if (num == 199999)
					break;
			}
		}
		file1.close();
		
		return veData;	//返回数据

	} while (false);

	return std::vector<double>();

}

void CDSP::DemDataFromFile()
{
	std::vector<double> dbData = GetDataFromFile();
	int iDataLen = dbData.size();
	if (iDataLen == 0)
		return;
	double *pData = new double[dbData.size()]();
	if (pData == NULL)
		return;
	for (int i = 0; i < iDataLen; i++)
	{
		pData[i] = dbData[i]/(1<<15);
	}

	QAM_Demod(pData, iDataLen, 100e6, 6.5e6, 1e6);

	delete[] pData;
}
