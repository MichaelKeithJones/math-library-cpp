#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

const long double PI = 3.141592653589793238L;

//--------------------------------------------------------------------------------------------------
// 3D Vector Class - Handles any vector of dimensions 2 or greater.
//--------------------------------------------------------------------------------------------------

template <typename T>
struct Vector
{
	int Dimension;
	T * values;

	//--------------------------------------------------------------------------------------------------
	// Constructor - A dimension of 2 or greater must be supplied. The default values of any matrix is 
	// an identity matrix.
	//--------------------------------------------------------------------------------------------------

	Vector(int dimension)
	{
		Dimension = dimension;
		values = new T[Dimension];

		for (int i = 0; i < Dimension; i++)
			values[i] = 0;

	}

	Vector(int dimension, const T & number)
	{
		Dimension = dimension;
		values = new T[Dimension];

		for (int i = 0; i < Dimension; i++)
			values[i] = number;
	}

	//--------------------------------------------------------------------------------------------------
	// Copy Constructor - Necessary when working with objects that contain pointers as member objects
	//--------------------------------------------------------------------------------------------------

	Vector(const Vector<T> & vector)
	{
		Dimension = vector.Dimension;
		values = new T[Dimension];

		for (int i = 0; i < Dimension; i++)
			values[i] = vector.values[i];
	}

	//--------------------------------------------------------------------------------------------------
	// Destructor - Frees dynamically allocated memory that "values" points too
	//--------------------------------------------------------------------------------------------------

	~Vector()
	{
		delete[] values; values = NULL;
	}

	//--------------------------------------------------------------------------------------------------
	// [] Operator Functions - Allows direct access to the data without having to call "values" member
	//--------------------------------------------------------------------------------------------------

	const T& operator[] (int index) const
	{
		return values[index];
	}

	T& operator[] (int index)
	{
		return values[index];
	}

	//--------------------------------------------------------------------------------------------------
	/* Vector Assignment:

		| 0  1  2  3 |	= | 0  1  2  3 |
	
	*/
	//--------------------------------------------------------------------------------------------------

	void operator= (const Vector<T> & vector)
	{
		Dimension = vector.Dimension;
		values = new T[Dimension];

		for (int i = 0; i < Dimension; i++)
			values[i] = vector.values[i];
	}

	//--------------------------------------------------------------------------------------------------
	/* Vector Addition:

		| 0  1  2  3 |	+ | 0  1  2  3 | = | 0+0  1+1  2+2  3+3 |

	*/
	//--------------------------------------------------------------------------------------------------

	Vector<T> operator+ (const Vector<T> & vector) const
	{
		Vector<T> temp;

		for (int i = 0; i < Dimension; i++)
			temp.values[i] = values[i] + vector.values[i];

		return temp;
	}

	//--------------------------------------------------------------------------------------------------
	/* Vector Subtraction:

		| 0  1  2  3 |	- | 0  1  2  3 | = | 0-0  1-1  2-2  3-3 |

	*/
	//--------------------------------------------------------------------------------------------------

	Vector<T> operator- (const Vector<T> & vector) const
	{
		Vector<T> temp;

		for (int i = 0; i < Dimension; i++)
			temp.values[i] = values[i] - vector.values[i];

		return temp;
	}

	//--------------------------------------------------------------------------------------------------
	/* Vector Multiplication (Scalar): 

		| 0  1  2  3 |	* A = | 0*A  1*A  2*A  3*A |

	*/
	//--------------------------------------------------------------------------------------------------

	Vector<T> operator* (const T & number) const
	{
		Vector<T> temp;

		for (int i = 0; i < Dimension; i++)
			temp.values[i] = values[i] * number;

		return temp;
	}
};

typedef Vector<float> vecf;

//--------------------------------------------------------------------------------------------------
// Matrix Class - Handles any square matrix of dimensions 2 or greater.
//--------------------------------------------------------------------------------------------------

template <typename T>
struct Matrix
{
	int Dimension;
	int Size;
	T * values;

	//--------------------------------------------------------------------------------------------------
	// Constructor - A dimension of 2 or greater must be supplied. The default values of any matrix is 
	// an identity matrix.
	//--------------------------------------------------------------------------------------------------

	Matrix(int dimension)
	{
		Dimension = dimension;
		Size = dimension * dimension;
		values = new T[Size];

		for (int i = 0; i < Size; i++)
		{
			if ((i % (Dimension + 1)) == 0)
				values[i] = 1;
			else
				values[i] = 0;
		}
	}

	//--------------------------------------------------------------------------------------------------
	// Copy Constructor - Necessary when working with objects that contain pointers as member objects
	//--------------------------------------------------------------------------------------------------

	Matrix(const Matrix & matrix)
	{
		Size = matrix.Size;
		Dimension = matrix.Dimension;
		values = new float[Size];

		for (int i = 0; i < Size; i++)
			values[i] = matrix.values[i];
	}

	//--------------------------------------------------------------------------------------------------
	// Destructor - Frees dynamically allocated memory that "values" points too
	//--------------------------------------------------------------------------------------------------

	~Matrix()
	{
		delete[] values; values = NULL;
	}
	
	//--------------------------------------------------------------------------------------------------
	// [] Operator Functions - Allows direct access to the data without having to call "values" member
	//--------------------------------------------------------------------------------------------------

	const T& operator[] (int index) const
	{
		return values[index];
	}

	T& operator[] (int index)
	{
		return values[index];
	}

	//--------------------------------------------------------------------------------------------------
	/* Matrix Assignment:

	| 0  4  8  12 |		| 0  4  8   12 |
	| 1  5  9  13 |	 =	| 1  5  9   13 |
	| 2  6  10 14 |		| 2  6  10  14 |
	| 3  7  11 15 |		| 3  7  11  15 |
	*/
	//--------------------------------------------------------------------------------------------------

	Matrix<T>& operator= (const Matrix<T> & matrix)
	{
		for (int i = 0; i < Size; i++)
			values[i] = matrix[i];

		return *this;
	}
	
	//--------------------------------------------------------------------------------------------------
	/* Matrix Addition:

	| 0  4  8  12 |		| 0  4  8   12 |    | 0 + 0  4 + 4  8 + 8	 12 + 12 |
	| 1  5  9  13 |	 +	| 1  5  9   13 | =	| 1 + 1  5 + 5  9 + 9    13 + 13 |
	| 2  6  10 14 |		| 2  6  10  14 |	| 2 + 2  6 + 6  10 + 10  14 + 14 |
	| 3  7  11 15 |		| 3  7  11  15 |	| 3 + 3  7 + 7  11 + 11  15 + 15 |
	*/
	//--------------------------------------------------------------------------------------------------

	Matrix<T> operator+ (const Matrix<T> & matrix)
	{
		Matrix<T> temp(Dimension);

		for (int i = 0; i < Size; i++)
			temp.values[i] = values[i] + matrix.values[i];

		return temp;
	}

	//--------------------------------------------------------------------------------------------------
	/* Matrix Subtraction:

	| 0  4  8  12 |		| 0  4  8   12 |    | 0 - 0  4 - 4  8 - 8	 12 - 12 |
	| 1  5  9  13 |	 -	| 1  5  9   13 | =	| 1 - 1  5 - 5  9 - 9    13 - 13 |
	| 2  6  10 14 |		| 2  6  10  14 |	| 2 - 2  6 - 6  10 - 10  14 - 14 |
	| 3  7  11 15 |		| 3  7  11  15 |	| 3 - 3  7 - 7  11 - 11  15 - 15 |
	*/
	//--------------------------------------------------------------------------------------------------

	Matrix<T> operator- (const Matrix<T> & matrix)
	{
		Matrix<T> temp(Dimension);

		for (int i = 0; i < Size; i++)
			temp.values[i] = values[i] - matrix.values[i];

		return temp;
	}

	//--------------------------------------------------------------------------------------------------
	/* Matrix Multiplication (Scalar):

	| 0  4  8  12 |			| 0 * A  4 * A  8 * A	12 * A |
	| 1  5  9  13 |	 * A =	| 1 * A  5 * A  9 * A   13 * A |
	| 2  6  10 14 |			| 2 * A  6 * A  10 * A  14 * A |
	| 3  7  11 15 |			| 3 * A  7 * A  11 * A  15 * A |
	*/
	//--------------------------------------------------------------------------------------------------

	Matrix<T> operator* (const T & number) const
	{
		Matrix<T> temp(Dimension);

		for (int i = 0; i < Size; i++)
			temp.values[i] = values[i] * number;

		return temp;
	}

	//--------------------------------------------------------------------------------------------------
	/* Matrix Multiplication (Matrices):

	| 0  4  8  12 |		| 0  4  8   12 |
	| 1  5  9  13 |	 *	| 1  5  9   13 | =	[i,j] = [i][j] + [i+4][j+1] +  [i+8][j+2] +  [i+12][j+3]
	| 2  6  10 14 |     | 2  6  10  14 |
	| 3  7  11 15 |		| 3  7  11  15 |
	-- j (0-4) , 16 , j+=4
	-- i (0-4) , 4 ,  i++
	*/
	//--------------------------------------------------------------------------------------------------

	Matrix<T> operator* (const Matrix<T> & input) const
	{
		Matrix<T> temp(Dimension);

		for (int i = 0; i < Size; i++)
			temp.values[i] = 0;

		for (int j = 0; j < Size; j += Dimension)
			for (int i = 0; i < Dimension; i++)
				for (int k = 0; k < Dimension; k++)
					temp.values[i + j] += values[i + k * Dimension] * input.values[j + k];

		return temp;
	}

	//--------------------------------------------------------------------------------------------------
	/* Matrix Transpose

	| 0  4  8  12 |		| 0   1   2   3 |
	| 1  5  9  13 |	 =	| 4   5   6   7 |
	| 2  6  10 14 |     | 8   9   10  11 |
	| 3  7  11 15 |		| 12  13  14  15 |

	*/
	//--------------------------------------------------------------------------------------------------

	Matrix<T> transpose() const
	{
		Matrix<T> temp(Dimension);

		for (int j = 0; j < Dimension; j++)
			for (int i = 0; i < Dimension; i++)
				temp.values[i * Dimension + j] = values[i + j * Dimension];

		return temp;
	}

	//--------------------------------------------------------------------------------------------------
	/* Matrix Determinant:

	| 0  2 | = (0 * 3) - (2 * 1)
	| 1  3 |


	| 0  3  6 |
	| 1  4  7 | = 0*det(4*8-7*5) - 3*det(1*8-7*2) + 6*det(1*5-4*2)
	| 2  5  8 |


	| 0  4  8  12 |
	| 1  5  9  13 | = 0*det(5,6,7,9,10,11,13,14,15) + 4*det(1,2,3,9,10,11,13,14,15) - 8*det(1,2,3,5,6,7,13,14,15) + 12*det(1,2,3,5,6,7,9,10,11)
	| 2  6  10 14 |
	| 3  7  11 15 |

		... 

	And so on recursively to other dimensions as long as the matrix is a square matrix of 2 or more in dimensions.
	"Determinant" function is a "kicker" function for the "det" sub-recursive function.

	*/
	//--------------------------------------------------------------------------------------------------

	T det(T * matrix, int size)
	{
		T result = 0;

		if (size == 2)
		{
			return (matrix[0] * matrix[3] - matrix[2] * matrix[1]);
		}

		float * submat = new float[(size - 1) * (size - 1)];
		int iterator = 0;
		int sign = 1;

		for (int j = 0; j < (size * size); j += size)
		{
			for (int i = 0; i < (size * size); i++)
			{
				if (i > j && i < (j + size)) {}
				else
					if (i % size != 0)
					{
						submat[iterator] = matrix[i];
						iterator++;
					}
			}
			iterator = 0;
			result += (matrix[j] * det(submat, size - 1)) * sign;
			sign *= -1;
		}

		delete[] submat;

		return result;
	}

	T determinant()
	{
		return det(values, Dimension);
	}

	//--------------------------------------------------------------------------------------------------
	/* Matrix Inverse (Row Major Method)

	| 0  4  8  12 |
	| 1  5  9  13 |
	| 2  6  10 14 |
	| 3  7  11 15 |

	Note: Row major method used here IAW: https://www.youtube.com/watch?v=nNOz6Ez8Fn4 but final answer saved in column major order for compatibility other functions and OpenGL

	*/
	//--------------------------------------------------------------------------------------------------

	Matrix<T> inverse()
	{
		Matrix<T> temp(Dimension);
		Matrix<T> submat(Dimension - 1);

		T determinant_inverse = determinant();
		determinant_inverse = T(1.0) / determinant_inverse;

		int sign = 1;
		int submat_iterator = 0;
		int temp_matrix_iterator = 0;

		for (int k = 0; k < Size; k += Dimension)
		{
			for (int j = 0; j < Dimension; j++)
			{
				for (int i = 0; i < Size; i++)
				{
					if (i >= k && i < (k + Dimension)) {}
					else
						if ((i - j) % Dimension != 0)
						{
							submat[submat_iterator] = values[i];
							submat_iterator++;
						}
				}
				temp[temp_matrix_iterator] = submat.determinant() * sign;
				cout << sign << " ";
				sign *= -1;


				submat_iterator = 0;
				temp_matrix_iterator++;
			}
			cout << endl;
			if((Dimension % 2) == 0) sign *= -1;
		}

		temp = temp.transpose();
		temp = temp * determinant_inverse;

		return temp;
	}
};

typedef Matrix<float> matf;

//--------------------------------------------------------------------------------------------------
// Print Vectors and Matrixes
//--------------------------------------------------------------------------------------------------

void print_vector(vecf temp)
{
	cout << setprecision(3);
	cout << "Vector: ";
	for (int i = 0; i < temp.Dimension; i++)
		cout << temp[i] << " ";
	cout << endl << endl;
}

void print_matrix(const matf & temp)
{
	cout << setprecision(3);
	cout << "Matrix: " << endl;

	for (int i = 0; i < temp.Dimension; i++)
	{
		for (int j = 0; j < (temp.Size); j+= temp.Dimension)
		{
			cout << temp.values[j+i] << " ";
		}
		cout << endl;
	}

	cout << endl;
}


//--------------------------------------------------------------------------------------------------
/*Projection Matrix Generator (Column Major Order):

	| 0  4  8  12 |
	| 1  5  9  13 |
	| 2  6  10 14 |
	| 3  7  11 15 |

*/
//--------------------------------------------------------------------------------------------------

template <typename T>
void genProjectionMatrix(T * temp, T FOV = 45.0f, T nearPlane = 0.1f, T farPlane = 100.0f, T imageWidth = 800.00, T imageHeight = 600.00)
{
	T fieldOfView;
	T aspectRatio;
	T top, bottom, left, right;

	fieldOfView = FOV * PI / 180;
	aspectRatio = imageWidth / imageHeight;

	top = tan(fieldOfView * 0.5) * nearPlane;
	right = top * aspectRatio;
	bottom = top * -1.0f;
	left = right * -1.0f;

	temp[0] = (2.0f * nearPlane) / (right - left);
	temp[5] = (2.0f * nearPlane) / (top - bottom);
	temp[8] = (right + left) / (right - left);
	temp[9] = (top + bottom) / (top - bottom);
	temp[10] = -(farPlane + nearPlane) / (farPlane - nearPlane);
	temp[11] = -1.0f;
	temp[14] = (-2.0 * farPlane * nearPlane) / (farPlane - nearPlane);
	temp[15] = 0.0f;
}

//--------------------------------------------------------------------------------------------------
// Cross Product:
//--------------------------------------------------------------------------------------------------

template <typename T>
T * cross(T *& v1, T *& v2)
{
	T * temp = new T[3];

	temp[0] = v1[1] * v2[2] - v1[2] * v2[1];
	temp[1] = v1[2] * v2[0] - v1[0] * v2[2];
	temp[2] = v1[0] * v2[1] - v1[1] * v2[0];

	return temp;
}

//--------------------------------------------------------------------------------------------------
// Dot Product:
//--------------------------------------------------------------------------------------------------

template <typename T>
T dot(T *& v1, T *& v2)
{
	T results = 0;

	for (int i = 0; i < 3; i++)
		results += (v1[i] * v2[i]);

	return results;
}

//--------------------------------------------------------------------------------------------------
// Length:
//--------------------------------------------------------------------------------------------------

template <typename T>
T length(T *& v1)
{
	T results = 0;

	for (int i = 0; i < 3; i++)
		results += (v1[i] * v1[i]);
	results = sqrt(results);

	return results;
}

//--------------------------------------------------------------------------------------------------
// Normalization:
//--------------------------------------------------------------------------------------------------

template <typename T>
void normalize(T *& v1)
{
	T len = length(v1);
	if (len > 0)
	{
		len = 1 / len;
		for (int i = 0; i < 3; i++)
			v1[i] *= len;
	}
}

//--------------------------------------------------------------------------------------------------
/* View Matrix Generator(Column Major Order) Right Hand Rule

From	-> To		-> Reference
Pos	-> Target	-> Up
Eye	-> Target	-> Up
Eye	-> Center	-> Up

*/
//--------------------------------------------------------------------------------------------------

template <typename T>
void genViewMatrix(T * temp, T * from, T * to, T * reference)
{
	T * forward, *right, *up;

	normalize(reference);

	forward = new T[3];
	for (int i = 0; i < 3; i++)
		forward[i] = from[i] - to[i];
	normalize(forward);

	right = cross(reference, forward);
	normalize(right);

	up = cross(forward, right);

	temp[0] = right[0];		temp[4] = right[1];		temp[8] = right[2];			temp[12] = -1 * dot(right, from);
	temp[1] = up[0];		temp[5] = up[1];		temp[9] = up[2];			temp[13] = -1 * dot(up, from);
	temp[2] = forward[0];	temp[6] = forward[1];	temp[10] = forward[2];		temp[14] = -1 * dot(forward, from);

	delete[] forward;
	delete[] right;
	delete[] up;
}

//--------------------------------------------------------------------------------------------------
// Main
//--------------------------------------------------------------------------------------------------

int main()
{
	float pos[3] = { 4, 3, -3 };
	float target[3] = {0, 0, 0};
	float up[3] = { 0, 1, 0 };


	cout << fixed;

	matf projection(4), view(4), model(4), mvp(4);
	genProjectionMatrix(projection.values);
		print_matrix(projection);

	genViewMatrix(view.values, pos, target, up);
		print_matrix(view);

	print_matrix(model);

	mvp = projection * view * model;

	print_matrix(mvp);

	cin.get();
	cin.get();

	return 0;
}