// This is the main DLL file.
#pragma once

#include "Forecast.h"

using namespace std;

extern "C" _declspec(dllexport) double* ComputeForecast(int forecastMethod, double* x, double* y, int inputLength, double* future, int futureLength){

	std::vector<double> xVector(x, x + inputLength);
	std::vector<double> yVector(y, y + inputLength);
	std::vector<double> futureVector(future, future + futureLength);

	DataArps dataArps(forecastMethod, xVector, yVector);

	std::vector<double> resultVector = dataArps.compute(futureVector);

	double* result = new double[futureLength];
	for (int i = 0; i < futureLength; i++){
		result[i] = resultVector[i];
	}
	return result;
}

extern "C" _declspec(dllexport) int ReleaseMemory(double* input){
	delete[] input;
	return 1;
}
