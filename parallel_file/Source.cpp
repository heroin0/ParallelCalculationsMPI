#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <mpi.h>

#include <ctime>
#ifdef _WIN32
#include <chrono>
#include <ppl.h>
using namespace std::chrono;
#endif

using namespace std;
class GravimetryTask
{
	private:
	vector<vector<vector<double> > >	firstSurfacePoints, secondSurfacePoints;
	vector<vector<double> > b,result;
	double dx;
	double dy;
	double gravityConst;// 6.67408e-11;
	int gridPointsAmount;
	int sqaredGridPointsAmount;
	double myEpsilon;

public:
	GravimetryTask(const string firstSurfaceFile,const string secondSurfaceFile,const string rightPartFile, const string resutlFile)
	{
		dx = 204.082;
		dy = 204.082;
		gravityConst = 6.6708e-3;// 6.67408e-11;
		gridPointsAmount = 50;
		sqaredGridPointsAmount = 2500;
		myEpsilon = 0.005;
		const int dimensions = 3;

		readDataFile(firstSurfaceFile, dimensions, firstSurfacePoints);
		readDataFile(secondSurfaceFile, dimensions, secondSurfacePoints);
		readDataFile(rightPartFile, dimensions, 2, b);
#ifdef _WIN32
		pplParallelEquationSolver(result);
#endif 
		printResultToFile(result, resutlFile);
	}
	//IO methods
	void readDataFile(string fileName, int dimensions, vector<vector<vector<double> > > &result)
	{
		vector<vector<vector<double> > >  tmpVector(gridPointsAmount, vector<vector<double> >(gridPointsAmount, vector<double>(dimensions)));
		ifstream  myFile;
		myFile.open(fileName.c_str());
		for (int i = 0; i < gridPointsAmount; i++)
		{
			for (int j = 0; j < gridPointsAmount; j++)
			{
				for (int k = 0; (k < dimensions) && (!myFile.eof()); k++)
				{
					myFile >> tmpVector[j][i][k];
				}
			}
		}
		myFile.close();
		result.swap(tmpVector);
	}

	void readDataFile(string fileName, int dimensions, int targetCollar, vector<vector<double> > &result)//should be faster
	{
		vector<vector<double> >  tmpVector(sqaredGridPointsAmount, vector<double>(1,0));
		ifstream  myFile;
		double dump;
		myFile.open(fileName.c_str());
		for (int i = 0; i < sqaredGridPointsAmount; i++)
		{
			for (int j = 0; j < dimensions&& (!myFile.eof()); j++)
			{
				if (j == targetCollar)
					myFile >> tmpVector[i][0];
				else
					myFile >> dump;
			}
		}
		myFile.close();
		swap(tmpVector,result);
		std::cout << "pew" << endl;
	}

	void printResultToFile(vector<vector<double> > z, string fileName)
	{
		ofstream myFile;
		myFile.open(fileName.c_str());
		for (int i = 0; i < gridPointsAmount; i++)
			for (int j = 0; j < gridPointsAmount; j++)
				myFile << firstSurfacePoints[i][j][0] << " " << firstSurfacePoints[i][j][1] << " " << z[j + i*gridPointsAmount][0] << "\n";
		myFile.close();
	}
	//Complicated math methods
#ifdef _WIN32
	void pplParallelEquationSolver(vector<vector<double> > &result)
	{
		vector<vector<double> > A, ATransposed, AonATransposed, ATranspozedOnB;
		vector<vector<double> > z(sqaredGridPointsAmount, vector<double>(1, 0)), tempMatrix1, tempMatrix2;
		double numerator, denominator;
		calculateMatrix(A);
		transpose(A, ATransposed);
		pplParallelMultiplicate(A, ATransposed, AonATransposed);
		pplParallelMultiplicate(ATransposed, b, ATranspozedOnB);
		double normB = vectorNorm(b);

		pplParallelMultiplicate(A, z, tempMatrix1);
		subtraction(tempMatrix1, b, tempMatrix2);
		double stoppingPoint = vectorNorm(tempMatrix2) / normB;


		int iterationNumber = 0;
		while (stoppingPoint >= myEpsilon)
		{
			iterationNumber++;
			pplParallelMultiplicate(AonATransposed, z, tempMatrix1);
			subtraction(tempMatrix1, ATranspozedOnB, tempMatrix2);
			numerator = pow(vectorNorm(tempMatrix2), 2);

			pplParallelMultiplicate(A, tempMatrix2, tempMatrix1);
			denominator = pow(vectorNorm(tempMatrix1), 2);

			pplParallelMultiplicate(A, z, tempMatrix1);
			subtraction(tempMatrix1, b, tempMatrix2);
			pplParallelMultiplicate(ATransposed, tempMatrix2, tempMatrix1);
			multiplicate(tempMatrix1, numerator / denominator, tempMatrix2);
			subtraction(z, tempMatrix2, tempMatrix1);
			z = tempMatrix1;

			pplParallelMultiplicate(A, z, tempMatrix1);
			subtraction(tempMatrix1, b, tempMatrix2);
			stoppingPoint = vectorNorm(tempMatrix2) / normB;
			if (iterationNumber % 100 == 0)
				cout << iterationNumber << ") " << stoppingPoint << endl;
		}
		cout << "Total iterations:" << iterationNumber << endl;
		swap(z, result);	
	}
#endif // _WIN32
	void calculateMatrix(vector<vector<double> > &result)
	{
		vector<vector<double> > tmpVector(gridPointsAmount*gridPointsAmount, vector<double>(gridPointsAmount*gridPointsAmount));

		for (int k = 0; k < gridPointsAmount; k++)
		{
			for (int l = 0; l < gridPointsAmount; l++)
			{
				for (int i = 0; i < gridPointsAmount; i++)//можно, конечно, использовать 1 брутальный цикл по 2500 строк.

					for (int j = 0; j < gridPointsAmount; j++)
					{
						tmpVector[gridPointsAmount*k + l][gridPointsAmount*i + j] = calculateIntegral(firstSurfacePoints[k][l][0], firstSurfacePoints[k][l][1], firstSurfacePoints[i][j][0], firstSurfacePoints[i][j][1], firstSurfacePoints[i][j][2], secondSurfacePoints[i][j][2]);
					}
			}
		}
		result.swap(tmpVector);
	}

	double calculateIntegral(double &xK, double &yL, double &xI, double &yJ, double &H1, double &H2)
	{
		double result = (1 / sqrt((xK - xI)*(xK - xI) + (yL - yJ)*(yL - yJ) + H1*H1));
		result -= (1 / sqrt((xK - xI)*(xK - xI) + (yL - yJ)*(yL - yJ) + H2*H2));
		result = gravityConst*dx*dy*result;
		return result;
	}

	//Matrix&vector arithmetics metods
	void transpose(vector<vector<double> > &inputMatrix, vector<vector<double> > &result)
	{
		const int m = inputMatrix.size(), n = inputMatrix[0].size();
		vector<vector<double> > tmpVector(m, vector<double>(n));
		//так себе транспозиция
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				tmpVector[i][j] = inputMatrix[j][i];
		swap(tmpVector, result);
	}

	void multiplicate(vector<vector<double> > &firstMatrix, vector<vector<double> > &secondMatrix, vector<vector<double> > &result)
	{
		const int m1 = firstMatrix.size(), m2 = secondMatrix.size(), n1 = firstMatrix[0].size(), n2 = secondMatrix[0].size();
		vector<vector<double> > tmpMatrix(m1, vector<double>(n2, 0));
		for (int i = 0; i < m1; i++)
			for (int j = 0; j < n2; j++)
				for (int k = 0; k < m2; k++)
					tmpMatrix[i][j] += firstMatrix[i][k] * secondMatrix[k][j];
		swap(tmpMatrix, result);
	}
	//only for use with Parallel Patterns Library
#ifdef _WIN32
	void pplParallelMultiplicate(vector<vector<double> > &firstMatrix, vector<vector<double> > &secondMatrix, vector<vector<double> > &result)
	{
		const int m1 = firstMatrix.size(), n1 = firstMatrix[0].size();
		const int m2 = secondMatrix.size(), n2 = secondMatrix[0].size();
		vector<vector<double> > tmpMatrix(m1, vector<double>(n2, 0));
		Concurrency::parallel_for(0, m1, [&](int i)
		{
			for (int j = 0; j < n2; j++)
				for (int k = 0; k < m2; k++)
					tmpMatrix[i][j] += firstMatrix[i][k] * secondMatrix[k][j];
		});
		swap(tmpMatrix, result);
	}
#endif
	void subtraction(vector<vector<double> > &firstMatrix, vector<vector<double> > &secondMatrix, vector<vector<double> > &result)
	{
		int m = firstMatrix.size(), n = firstMatrix[0].size();
		vector<vector<double> > tmpMatrix(m, vector<double>(n, 0));
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < n; j++)
				tmpMatrix[i][j] = firstMatrix[i][j] - secondMatrix[i][j];
		}
		swap(tmpMatrix, result);
	}

	double vectorNorm(vector<vector<double> > &matrix)
	{
		double result = 0;
		for (unsigned int i = 0; i < matrix.size(); i++)
			result += matrix[i][0] * matrix[i][0];
		return sqrt(result);
	}

	void multiplicate(vector<vector<double> > &matrix, double number, vector<vector<double> > &result)
	{
		const int m = matrix.size(), n = matrix[0].size();
		vector<vector<double> > tmpVector(m, vector<double>(n));
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				tmpVector[i][j] = matrix[i][j] * number;
		swap(result, tmpVector);
	}
};

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	steady_clock::time_point t1 = steady_clock::now();
	GravimetryTask g("hh1.dat", "hh2.dat", "f_1sq.dat","result.dat");
	steady_clock::time_point t2 = steady_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	cout << "It took me " << time_span.count() << " seconds." << endl;;
	MPI_Finalize();
	return 0;
}