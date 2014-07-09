// This is the main DLL file.
#pragma once

#include "Forecast.h"

using namespace std;

extern "C" _declspec(dllexport) double* ComputeArpsForecast(int forecastMethod, double* x, double* y, int inputLength, double* future, int futureLength){

	std::vector<double> xVector(x, x + inputLength);
	std::vector<double> yVector(y, y + inputLength);
	std::vector<double> futureVector(future, future + futureLength);

	DataArps dataArps(forecastMethod, xVector, yVector);

	std::vector<double> resultVector = dataArps.computeForecast(futureVector);

	double* result = new double[futureLength];
	for (int i = 0; i < futureLength; i++){
		result[i] = resultVector[i];
	}
	return result;
}

extern "C" _declspec(dllexport) double ComputeArpsEur(int forecastMethod, double* x, double* y, int inputLength, double Qi, double qf){

	std::vector<double> xVector(x, x + inputLength);
	std::vector<double> yVector(y, y + inputLength);

	DataArps dataArps(forecastMethod, xVector, yVector);

	return dataArps.computeEur(Qi, qf);
}

extern "C" _declspec(dllexport) double* ComputeDuongForecast(int forecastMethod, double* x, double* y, int inputLength, double* future, int futureLength){

	std::vector<double> xVector(x, x + inputLength);
	std::vector<double> yVector(y, y + inputLength);
	std::vector<double> futureVector(future, future + futureLength);

	DataDuong dataDuong(forecastMethod, xVector, yVector);

	std::vector<double> resultVector = dataDuong.computeForecast(futureVector);

	double* result = new double[futureLength];
	for (int i = 0; i < futureLength; i++){
		result[i] = resultVector[i];
	}
	return result;
}

extern "C" _declspec(dllexport) double ComputeDuongEur(int forecastMethod, double* x, double* y, int inputLength, double tf){

	std::vector<double> xVector(x, x + inputLength);
	std::vector<double> yVector(y, y + inputLength);

	DataDuong dataDuong(forecastMethod, xVector, yVector);

	return dataDuong.computeEur(tf);
}

extern "C" _declspec(dllexport) double* ComputeSEPDForecast(double* x, double* y, int inputLength, double* future, int futureLength){

	std::vector<double> xVector(x, x + inputLength);
	std::vector<double> yVector(y, y + inputLength);
	std::vector<double> futureVector(future, future + futureLength);

	DataSEPD dataSEPD(xVector, yVector);

	std::vector<double> resultVector = dataSEPD.computeForecast(futureVector);

	double* result = new double[futureLength];
	for (int i = 0; i < futureLength; i++){
		result[i] = resultVector[i];
	}
	return result;
}

extern "C" _declspec(dllexport) double ComputeSEPDEur(double* x, double* y, int inputLength){

	std::vector<double> xVector(x, x + inputLength);
	std::vector<double> yVector(y, y + inputLength);

	DataSEPD dataSEPD(xVector, yVector);

	return dataSEPD.ComputeEur();
}

extern "C" _declspec(dllexport) int ReleaseMemory(double* input){
	delete[] input;
	return 1;
}
