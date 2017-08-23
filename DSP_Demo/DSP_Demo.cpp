// DSP_Demo.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include "DSP.h"
#include "fir.h"
#include <time.h>
#include <fstream>

int _tmain(int argc, _TCHAR* argv[])
{
	std::vector<double> pOut;
	double xn[3] = {3.0, 4.0, 5.0};
	double hn[3] = { 1.0, 2.0, 3.0};
	/*double xn[3] = {0.0, 4.0, 3.0};
	double hn[2] = { 2.0, 1.0 };*/
	CDSP Obj_dsp;
	pOut = Obj_dsp.FastConv(xn, 3, hn, 3);
	for (int i = 0; i < 5; i++)
		printf("%.2f\t", pOut[i]);
	//Obj_dsp.QAM_Gen(100e6, 1e6, 400);
	Obj_dsp.DemDataFromFile();

	enWindType  wType = WINDOW_HAMMING;
	/*enFiltType  fType = FILTER_LOWPASS;
	double  fs = 100e6,
			fpass = 1e6,
			apass = -3,
			fstop = 2e6,
			astop = -20;

	CFIR fir(fType, wType);*/
	enFiltType fType = FILTER_BANDSTOP;
		double	fs = 1000,
				fpass1 = 100,
				apass1 = -3,
				fstop1 = 200,
				fstop2 = 300,
				astop1 = -20,
				fpass2 = 400,
				apass2 = -3;
	CFIR fir(fType, wType);
	fir.setParams(fs, fpass1, apass1, fstop1, fstop2, astop1, fpass2, apass2);
	fir.design();
	//fir.dispInfo();
	std::cout << std::endl;
	/*std::vector<double>	coefs = fir.getCoefs();
	std::cout << coefs.size() << std::endl;*/



	system("pause");
	return 0;
}

