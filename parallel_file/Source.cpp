#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <mpi.h>
#include <chrono>
#include <ctime>
#include <ppl.h>
using namespace std;
//метод наискорейшего спуска
using namespace std::chrono;

class GravimetryTask
{
private:
	vector<vector<vector<double> > >	firstSurfacePoints, secondSurfacePoints, rightPartPoints;

	double dx;
	double dy;
	double gravityConst;// 6.67408e-11;
	int gridPointsAmount;
	int sqaredGridPointsAmount;
public:
	GravimetryTask(string firstSurfaceFile, string secondSurfaceFile, string rightPartFile)
	{
		dx = 204.082;
		dy = 204.082;
		gravityConst = 6.6708e-3;// 6.67408e-11;
		gridPointsAmount = 50;
		sqaredGridPointsAmount = 2500;//константные переменные класса
		const int dimensions = 3;

		readDataFile(firstSurfaceFile, dimensions, firstSurfacePoints);
		readDataFile(secondSurfaceFile, dimensions, secondSurfacePoints);
		readDataFile(secondSurfaceFile, dimensions, rightPartPoints);
		vector<vector<double> > equationMatrix, transposedMatrix, AonATransposed, z(1, vector<double>(sqaredGridPointsAmount, 0));
		calculateMatrix(equationMatrix);
		transpose(equationMatrix, transposedMatrix);
		parallelMultiplicate(equationMatrix, transposedMatrix, AonATransposed);
		//multiplicate(equationMatrix, transposedMatrix, AonATransposed);
	}

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

	void parallelMultiplicate(vector<vector<double> > &firstMatrix, vector<vector<double> > &secondMatrix, vector<vector<double> > &result)
	{
		const int m1 = firstMatrix.size(), m2 = secondMatrix.size(), n1 = firstMatrix[0].size(), n2 = secondMatrix[0].size();
		vector<vector<double> > tmpMatrix(m1, vector<double>(n2, 0));
		Concurrency::parallel_for(0, m1, [&](int i)
		{
			for (int j = 0; j < n2; j++)
				for (int k = 0; k < m2; k++)
					tmpMatrix[i][j] += firstMatrix[i][k] * secondMatrix[k][j];
		});
		swap(tmpMatrix, result);
	}

	void subtraction(vector<vector<double> > &firstMatrix, vector<vector<double> > &secondMatrix, vector<vector<double> > &result)
	{
		int m = firstMatrix.size(), n = firstMatrix[0].size();
		vector<vector<double> > tmpMatrix(m, vector<double>(n, 0));
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				tmpMatrix[i][j] = firstMatrix[i][j] - secondMatrix[i][j];
		swap(tmpMatrix, result);
	}

	double vectorNorm(vector<vector<double> > &matrix)
	{
		double result = 0;
		for (unsigned int i = 0; i < matrix.size(); i++)
			result += matrix[0][i] * matrix[0][i];
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

	void readDataFile(string fileName, int dimensions, vector<vector<vector<double> > > &result)//should be faster
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

	void calculateMatrix(vector<vector<double> > &result)
	{
		vector<vector<double> > tmpVector(gridPointsAmount*gridPointsAmount, vector<double>(gridPointsAmount*gridPointsAmount));

		for (int k = 0; k < gridPointsAmount; k++)
		{
			for (int l = 0; l < gridPointsAmount; l++)
			{
				for (int i = 0; i < gridPointsAmount; i++)//можно, конечно, использовать 1 брутальный цикл по 2500 строк.

					for (int j = 0; j < gridPointsAmount; j++)
					{//gromozdko 
						tmpVector[gridPointsAmount*k + l][gridPointsAmount*i + j] = calculateIntegral(firstSurfacePoints[k][l][0], firstSurfacePoints[k][l][1], firstSurfacePoints[i][j][0], firstSurfacePoints[i][j][1], firstSurfacePoints[i][j][2], secondSurfacePoints[i][j][2]);
					}
			}
		}

		result.swap(tmpVector);
	}

	double calculateIntegral(double &xK, double &yL, double &xI, double &yJ, double &H1, double &H2)
	{
		double result = (1 / sqrt((xK - xI)*(xK - xI) + (yL - yJ)*(yL - yJ) + H1*H1));//do not use pow cause of loose precision
		result -= (1 / sqrt((xK - xI)*(xK - xI) + (yL - yJ)*(yL - yJ) + H2*H2));
		result = gravityConst*dx*dy*result;
		return result;
	}
};

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	steady_clock::time_point t1 = steady_clock::now();
	GravimetryTask g("hh1.dat", "hh2.dat", "f_1sq.dat");
	steady_clock::time_point t2 = steady_clock::now();
	duration<double> time_span = duration_cast<duration<double> >(t2 - t1);
	cout << "It took me " << time_span.count() << " seconds." << endl;;
	MPI_Finalize();
	return 0;
}

//GravimetryTask g("f_1sq.dat","f_2sq_cl.dat"); //(Стас) это правые части, use hh files instead, Luke!