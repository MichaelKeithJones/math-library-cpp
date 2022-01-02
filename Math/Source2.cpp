#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

const long double PI = 3.141592653589793238L;

//--------------------------------------------------------------------------------------------------
// 4D Matrix Structure
//--------------------------------------------------------------------------------------------------

struct Matrix
{
	int Size;
	int Dimension;
	float * values;

	Matrix()
	{
		Size = 16;
		Dimension = 4;
		values = new float[Size];

		for (int i = 0; i < Size; i++)
		{
			if ((i % (Dimension + 1)) == 0)
				values[i] = 1;
			else
				values[i] = 0;
		}
	}

	Matrix(const Matrix & matrix)
	{
		Size = matrix.Size;
		Dimension = matrix.Dimension;
		values = new float[Size];

		for (int i = 0; i < Size; i++)
			values[i] = matrix.values[i];
	}

	~Matrix()
	{
		delete[] values; values = NULL;
	}

	const float& operator[] (int index) const
	{
		return values[index];
	}

	float& operator[] (int index)
	{
		return values[index];
	}

	void operator= (const Matrix & matrix)
	{
		for (int i = 0; i < Size; i++)
			values[i] = matrix.values[i];
	}
};


Matrix operator+ (const Matrix & m1, const Matrix & m2)
{
	Matrix temp;

	for (int i = 0; i < m1.Size; i++)
		temp.values[i] = m1.values[i] + m2.values[i];

	return temp;
}

//--------------------------------------------------------------------------------------------------
// Print  Matrixes
//--------------------------------------------------------------------------------------------------

void print_matrix(Matrix temp)
{
	cout << setprecision(3);
	cout << "Matrix: " << temp[0] << " " << temp[4] << " " << temp[8] << " " << temp[12] << endl
		<< "        " << temp[1] << " " << temp[5] << " " << temp[9] << " " << temp[13] << endl
		<< "        " << temp[2] << " " << temp[6] << " " << temp[10] << " " << temp[14] << endl
		<< "        " << temp[3] << " " << temp[7] << " " << temp[11] << " " << temp[15] << endl;

	cout << endl;
}

int main()
{
	Matrix m;
	Matrix n;

	print_matrix(m);
	print_matrix(n);

	m = m + n;

	print_matrix(m);

	cin.get();
	cin.get();

	return 0;
}
