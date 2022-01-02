#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

const long double PI = 3.141592653589793238L;

//--------------------------------------------------------------------------------------------------
// 3D Vector Structure
//--------------------------------------------------------------------------------------------------

template <typename T>
struct vector3
{
	vector3();
	vector3(const T & number);
	vector3(T x, T y, T z);

	void operator= (const vector3<T> &);

	const T& operator[] (int index) const;
	T& operator[] (int index);

	vector3<T> operator+ (const vector3<T> & vector) const;
	vector3<T> operator- (const vector3<T> & vector) const;
	vector3<T> operator* (const T & number) const;

	vector3<T> cross(const vector3<T> & vector) const;
	T dot(const vector3<T> & vector) const;

	T length() const;
	void normalize();
	
private:
	T values[3] = {0,0,0};
};
typedef vector3<float> vec3f;

//--------------------------------------------------------------------------------------------------
// 4D Matrix Structure
//--------------------------------------------------------------------------------------------------

template <typename T>
struct matrix4
{
	void operator= (const matrix4<T> & matrix);

	const T& operator[] (int index) const;
	T& operator[] (int index);

	matrix4<T> operator+ (const matrix4<T> & matrix) const;
	matrix4<T> operator- (const matrix4<T> & matrix) const;
	matrix4<T> operator* (const matrix4<T> & input) const;
	matrix4<T> operator* (const T & number) const;

	matrix4<T> transpose() const;
	matrix4<T> inverse() const;

private:
	T values[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
};
typedef matrix4<float> mat4f;

//--------------------------------------------------------------------------------------------------
// Matrix Functions
//--------------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------------
/* 2D Matrix Determinant:

| 0  2 | = (0 * 3) - (2 * 1)
| 1  3 |

*/
//--------------------------------------------------------------------------------------------------

template <typename T>
T determinant2x2(const T & v0, const T & v1, const T & v2, const T & v3)
{
	return (v0*v3 - v2*v1);
}

//--------------------------------------------------------------------------------------------------
/* 3D Matrix Determinant:

| 0  3  6 |
| 1  4  7 | = 0*det(4*8-7*5) - 3*det(1*8-7*2) + 6*det(1*5-4*2)
| 2  5  8 |

*/
//--------------------------------------------------------------------------------------------------

template <typename T>
T determinant3x3(const T & v0, const T & v1, const T & v2, const T & v3, const T & v4, const T & v5, const T & v6, const T & v7, const T & v8)
{
	return (v0 * determinant2x2(v4,v5,v7,v8) - v3 * determinant2x2(v1,v2,v7,v8) + v6 * determinant2x2(v1,v2,v4,v5));
}

//--------------------------------------------------------------------------------------------------
/* 4D Matrix Determinant:

| 0  4  8  12 |
| 1  5  9  13 | = 0*det(5,6,7,9,10,11,13,14,15) + 4*det(1,2,3,9,10,11,13,14,15) - 8*det(1,2,3,5,6,7,13,14,15) + 12*det(1,2,3,5,6,7,9,10,11)
| 2  6  10 14 |
| 3  7  11 15 |

*/
//--------------------------------------------------------------------------------------------------

template <typename T>
T determinant4x4(const T & v0, const T & v1, const T & v2, const T & v3, const T & v4, const T & v5, const T & v6, const T & v7, const T & v8, const T & v9, const T & v10, const T & v11, const T & v12, const T & v13, const T & v14, const T & v15)
{
	return (v0 * determinant3x3(v5, v6, v7, v9, v10, v11, v13, v14, v15) + v4 * determinant3x3(v1, v2, v3, v9, v10, v11, v13, v14, v15) - v8 * determinant3x3(v1, v2, v3, v5, v6, v7, v13, v14, v15) + v12 * determinant3x3(v1, v2, v3, v5, v6, v7, v9, v10, v11));
}

//--------------------------------------------------------------------------------------------------
//Projection Matrix Generator (Column Major Order)
//--------------------------------------------------------------------------------------------------

void genProjectionMatrix(mat4f & temp, float FOV = 45.0f, float nearPlane = 0.1f, float farPlane = 100.0f, unsigned int imageWidth = 800, unsigned int imageHeight = 600)
{
	float fieldOfView;
	float aspectRatio;
	float top, bottom, left, right;
		
	fieldOfView = FOV * PI / 180;
	aspectRatio = (float)imageWidth / (float)imageHeight;
		
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
/* View Matrix Generator(Column Major Order) Right Hand Rule

 From	-> To		-> Reference
 Pos	-> Target	-> Up
 Eye	-> Target	-> Up
 Eye	-> Center	-> Up

*/
//--------------------------------------------------------------------------------------------------

void genViewMatrix(mat4f & temp, const vec3f & from, const vec3f & to, vec3f reference)
{
	vec3f right, forward, up;

	reference.normalize();

	forward = from - to;
	forward.normalize();

	right = reference.cross(forward);
	right.normalize();

	up = forward.cross(right);

	temp[0] = right[0];		temp[4] = right[1];		temp[8] = right[2];			temp[12] = -right.dot(from);
	temp[1] = up[0];		temp[5] = up[1];		temp[9] = up[2];			temp[13] = -up.dot(from);
	temp[2] = forward[0];	temp[6] = forward[1];	temp[10] = forward[2];		temp[14] = -forward.dot(from);
}

//--------------------------------------------------------------------------------------------------
// Print Vectors and Matrixes
//--------------------------------------------------------------------------------------------------

void print_vector(vec3f temp)
{
	cout << setprecision(3);
	cout << "Vector: " << temp[0] << " " << temp[1] << " " << temp[2] << " " << temp[3] << endl;

	cout << endl;
}

void print_matrix(mat4f temp)
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
	vec3f pos(4,3,-3), target(0,0,0), up(0,1,0);
	mat4f projection, view, model, mvp;

	genProjectionMatrix(projection);
	print_matrix(projection);

	genViewMatrix(view, pos, target, up);
	print_matrix(view);

	print_matrix(model);

	mvp = projection * view * model;

	print_matrix(mvp);

	mat4f temp;
	for (int i = 0; i < 16; i++) temp[i] = i;
	print_matrix(temp);

	temp = temp.transpose();
	print_matrix(temp);

	cin.get();
	cin.get();

	return 0;
}

//--------------------------------------------------------------------------------------------------
// 3D Vector Definitions
//--------------------------------------------------------------------------------------------------

template <typename T>
vector3<T>::vector3() {}

template <typename T>
vector3<T>::vector3(const T & number)
{
	values[0] = values[1] = values[2] = number;
}

template <typename T>
vector3<T>::vector3(T x, T y, T z)
{
	values[0] = x;
	values[1] = y;
	values[2] = z;
}

template <typename T>
void vector3<T>::operator= (const vector3<T> & vector)
{
	values[0] = vector.values[0];
	values[1] = vector.values[1];
	values[2] = vector.values[2];
}

template <typename T>
const T& vector3<T>::operator[] (int index) const
{ 
	return values[index];
}

template <typename T>
T& vector3<T>::operator[] (int index)
{ 
	return values[index];
}

template <typename T>
vector3<T> vector3<T>::operator+ (const vector3<T> & vector) const
{
	vector3<T> temp;

	temp.values[0] = values[0] + vector.values[0];
	temp.values[1] = values[1] + vector.values[1];
	temp.values[2] = values[2] + vector.values[2];

	return temp;
}

template <typename T>
vector3<T> vector3<T>::operator- (const vector3<T> & vector) const
{
	vector3<T> temp;

	temp.values[0] = values[0] - vector.values[0];
	temp.values[1] = values[1] - vector.values[1];
	temp.values[2] = values[2] - vector.values[2];

	return temp;
}

template <typename T>
vector3<T> vector3<T>::operator* (const T & number) const
{
	vector3<T> temp;

	temp.values[0] = values[0] * number;
	temp.values[1] = values[1] * number;
	temp.values[2] = values[2] * number;

	return temp;
}

template <typename T>
vector3<T> vector3<T>::cross(const vector3<T> & vector) const
{
	vector3<T> temp;

	temp.values[0] = values[1] * vector.values[2] - values[2] * vector.values[1];
	temp.values[1] = values[2] * vector.values[0] - values[0] * vector.values[2];
	temp.values[2] = values[0] * vector.values[1] - values[1] * vector.values[0];

	return temp;
}

template <typename T>
T vector3<T>::dot(const vector3<T> & vector) const
{ 
	return(values[0] * vector.values[0] + values[1] * vector.values[1] + values[2] * vector.values[2]); 
}

template <typename T>
T vector3<T>::length() const
{ 
	return sqrt(values[0] * values[0] + values[1] * values[1] + values[2] * values[2]); 
}

template <typename T>
void vector3<T>::normalize()
{
	T len = length();
	if (len > 0)
	{
		len = 1 / len;
		values[0] *= len;
		values[1] *= len;
		values[2] *= len;
	}
}

//--------------------------------------------------------------------------------------------------
/* 4D Matrix Definitions:

	Matrices are in column-major orientation in order to be inline with OpenGL:

	| 0  4  8  12 |
	| 1  5  9  13 |
	| 2  6  10 14 |
	| 3  7  11 15 |
*/
//--------------------------------------------------------------------------------------------------


template <typename T>
void matrix4<T>::operator= (const matrix4<T> & matrix)
{
	for (int i = 0; i < 16; i++)
		values[i] = matrix.values[i];
}

template <typename T>
const T& matrix4<T>::operator[] (int index) const
{ 
	return values[index]; 
}

template <typename T>
T& matrix4<T>::operator[] (int index)
{ 
	return values[index]; 
}

//--------------------------------------------------------------------------------------------------
/* 4D Matrix Addition:

	| 0  4  8  12 |		| 0  4  8   12 |    | 0 + 0  4 + 4  8 + 8	 12 + 12 |
	| 1  5  9  13 |	 +	| 1  5  9   13 | =	| 1 + 1  5 + 5  9 + 9    13 + 13 | 
	| 2  6  10 14 |		| 2  6  10  14 |	| 2 + 2  6 + 6  10 + 10  14 + 14 |
	| 3  7  11 15 |		| 3  7  11  15 |	| 3 + 3  7 + 7  11 + 11  15 + 15 |
*/
//--------------------------------------------------------------------------------------------------

template <typename T>
matrix4<T> matrix4<T>::operator+ (const matrix4<T> & matrix) const
{
	matrix4<T> temp;

	for (int i = 0; i < 16; i++)
		temp.values[i] = values[i] + matrix.values[i];

	return temp;
}

//--------------------------------------------------------------------------------------------------
/* 4D Matrix Subtraction:

| 0  4  8  12 |		| 0  4  8   12 |    | 0 - 0  4 - 4  8 - 8	 12 - 12 |
| 1  5  9  13 |	 -	| 1  5  9   13 | =	| 1 - 1  5 - 5  9 - 9    13 - 13 |
| 2  6  10 14 |		| 2  6  10  14 |	| 2 - 2  6 - 6  10 - 10  14 - 14 |
| 3  7  11 15 |		| 3  7  11  15 |	| 3 - 3  7 - 7  11 - 11  15 - 15 |
*/
//--------------------------------------------------------------------------------------------------

template <typename T>
matrix4<T> matrix4<T>::operator- (const matrix4<T> & matrix) const
{
	matrix4<T> temp;

	for (int i = 0; i < 16; i++)
		temp.values[i] = values[i] - matrix.values[i];

	return temp;
}

//--------------------------------------------------------------------------------------------------
/* 4D Matrix Multiplication (Scalar):

| 0  4  8  12 |			| 0 * A  4 * A  8 * A	12 * A |
| 1  5  9  13 |	 * A =	| 1 * A  5 * A  9 * A   13 * A |
| 2  6  10 14 |			| 2 * A  6 * A  10 * A  14 * A |
| 3  7  11 15 |			| 3 * A  7 * A  11 * A  15 * A |
*/
//--------------------------------------------------------------------------------------------------

template <typename T>
matrix4<T> matrix4<T>::operator* (const T & number) const
{
	matrix4<T> temp;

	for (int i = 0; i < 16; i++)
		temp.values[i] = values[i] * number;

	return temp;
}

//--------------------------------------------------------------------------------------------------
/* 4D Matrix Multiplication (Matrices):

| 0  4  8  12 |		| 0  4  8   12 |   
| 1  5  9  13 |	 *	| 1  5  9   13 | =	[i,j] = [i][j] + [i+4][j+1] +  [i+8][j+2] +  [i+12][j+3]
| 2  6  10 14 |     | 2  6  10  14 |
| 3  7  11 15 |		| 3  7  11  15 |	
										-- j (0-4) , 16 , j+=4
										-- i (0-4) , 4 ,  i++
*/
//--------------------------------------------------------------------------------------------------

template <typename T>
matrix4<T> matrix4<T>::operator* (const matrix4<T> & input) const
{
	matrix4<T> temp;

	for (int j = 0; j < 16; j+=4)
		for (int i = 0; i < 4; i++)
			temp.values[i+j] = (values[i] * input.values[j]) + (values[i+4] * input.values[j+1]) + (values[i+8] * input.values[j+2]) + (values[i+12] * input.values[j+3]);

	return temp;
}

//--------------------------------------------------------------------------------------------------
/* 4D Matrix Transpose

| 0  4  8  12 |		| 0   1   2   3 |
| 1  5  9  13 |	 =	| 4   5   6   7 |
| 2  6  10 14 |     | 8   9   10  11 |
| 3  7  11 15 |		| 12  13  14  15 |

*/
//--------------------------------------------------------------------------------------------------

template <typename T>
matrix4<T> matrix4<T>::transpose() const
{
	matrix4<T> temp;

	temp.values[0] = values[0];		temp.values[4] = values[1];		temp.values[8] = values[2];		temp.values[12] = values[3];
	temp.values[1] = values[4];		temp.values[5] = values[5];		temp.values[9] = values[6];		temp.values[13] = values[7];
	temp.values[2] = values[8];		temp.values[6] = values[9];		temp.values[10] = values[10];	temp.values[14] = values[11];
	temp.values[3] = values[12];	temp.values[7] = values[13];	temp.values[11] = values[14];	temp.values[15] = values[15];
	 
	return temp;
}

//--------------------------------------------------------------------------------------------------
/* 4D Matrix Inverse (Row Major Method)

| 0  4  8  12 |		
| 1  5  9  13 |	
| 2  6  10 14 |   
| 3  7  11 15 |

Note: Row major method used here IAW: https://www.youtube.com/watch?v=nNOz6Ez8Fn4 but final answer saved in column major order for compatibility other functions and OpenGL

*/
//--------------------------------------------------------------------------------------------------

template <typename T>
matrix4<T> matrix4<T>::inverse() const
{
	matrix4<T> temp;

	T det = determinant4x4(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], values[10], values[11], values[12], values[13], values[14], values[15]);
	det = T(1.0) / det;

	temp[0] = determinant3x3(values[5], values[6], values[7], values[9], values[10], values[11], values[13], values[14], values[15]);
	temp[4] = -determinant3x3(values[1], values[2], values[3], values[9], values[10], values[11], values[13], values[14], values[15]);
	temp[8] = determinant3x3(values[1], values[2], values[3], values[5], values[6], values[7], values[13], values[14], values[15]);
	temp[12] = -determinant3x3(values[1], values[2], values[3], values[5], values[6], values[7], values[9], values[10], values[11]);

	temp[1] = -determinant3x3(values[4], values[6], values[7], values[8], values[10], values[11], values[12], values[14], values[15]);
	temp[5] = determinant3x3(values[0], values[2], values[3], values[8], values[10], values[11], values[12], values[14], values[15]);
	temp[9] = -determinant3x3(values[0], values[2], values[3], values[4], values[6], values[7], values[12], values[14], values[15]);
	temp[13] = determinant3x3(values[0], values[2], values[3], values[4], values[6], values[7], values[8], values[10], values[11]);

	temp[2] = determinant3x3(values[4], values[5], values[7], values[8], values[9], values[11], values[12], values[13], values[15]);
	temp[6] = -determinant3x3(values[0], values[1], values[3], values[8], values[9], values[11], values[12], values[13], values[15]);
	temp[10] = determinant3x3(values[0], values[1], values[3], values[4], values[5], values[7], values[12], values[13], values[15]);
	temp[14] = -determinant3x3(values[0], values[1], values[3], values[4], values[5], values[7], values[8], values[9], values[11]);

	temp[3] = -determinant3x3(values[4], values[5], values[6], values[8], values[9], values[10], values[12], values[13], values[14]);
	temp[7] = determinant3x3(values[0], values[1], values[2], values[8], values[9], values[10], values[12], values[13], values[14]);
	temp[11] = -determinant3x3(values[0], values[1], values[2], values[4], values[5], values[6], values[12], values[13], values[14]);
	temp[15] = determinant3x3(values[0], values[1], values[2], values[4], values[5], values[6], values[8], values[9], values[10]);

	temp = temp.transpose();
	temp = temp * det;

	return temp;
}