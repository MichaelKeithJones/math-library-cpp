
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

const long double PI = 3.141592653589793238L;

//--------------------------------------------------------------------------------------------------
// 2D Vector Structure
//--------------------------------------------------------------------------------------------------

template <typename T>
struct vector2
{
	vector2() { values[0] = values[1] = 0; }
	vector2(const T & number) { values[0] = values[1] = number; }
	vector2(T x, T y);

	void operator= (const vector2<T>  &);

	const T& operator[] (int index) const { return values[index]; }
	T& operator[] (int index) { return values[index]; }

	vector2<T> operator+ (const vector2<T>  & vector) const;
	vector2<T> operator- (const vector2<T>  & vector) const;

	T length() const { return sqrt(values[0] * values[0] + values[1] * values[1]); }
	void normalize();

	T dot(const vector2<T>  & vector) const { return(values[0] * vector.values[0] + values[1] * vector.values[1]); }
	vector2<T> operator* (const T & number) const;
	T operator* (const vector2<T>  & vector) const;

private:
	T values[2];
};
typedef vector2<float> vec2f;

//--------------------------------------------------------------------------------------------------
// 3D Vector Structure
//--------------------------------------------------------------------------------------------------

template <typename T>
struct vector3
{
	vector3() { values[0] = values[1] = values[2] = 0;}
	vector3(const T & number) { values[0] = values[1] = values[2] = number;}
	vector3(T x, T y, T z);

	void operator= (const vector3<T> &);

	const T& operator[] (int index) const { return values[index]; }
	T& operator[] (int index) { return values[index]; }

	vector3<T> operator+ (const vector3<T> & vector) const;
	vector3<T> operator- (const vector3<T> & vector) const;

	T length() const { return sqrt(values[0] * values[0] + values[1] * values[1] + values[2] * values[2]); }
	void normalize();

	T dot(const vector3<T> & vector) const { return(values[0] * vector.values[0] + values[1] * vector.values[1] + values[2] * vector.values[2]); }
	vector3<T> operator* (const T & number) const;
	vector3<T> operator* (const vector3<T> & vector) const;

private:
	T values[3];
};
typedef vector3<float> vec3f;


//--------------------------------------------------------------------------------------------------
// 4D Vector Structure
//--------------------------------------------------------------------------------------------------

template <typename T>
struct vector4
{
	vector4()						{	values[0] = values[1] = values[2] = 0; values[3] = 1.0f; }
	vector4(const T & number)		{	values[0] = values[1] = values[2] = number; values[3] = 1.0f; }
	vector4(T x, T y, T z);

	void operator= (const vector4<T> &);

	const T& operator[] (int index) const	{	return values[index];	}
	T& operator[] (int index)	{	return values[index];	}

	vector4<T> operator+ (const vector4<T> & vector) const;
	vector4<T> operator- (const vector4<T> & vector) const;

	T length() const {	return sqrt(values[0] * values[0] + values[1] * values[1] + values[2] * values[2]);	}
	void normalize();

	T dot (const vector4<T> & vector) const	{	return(values[0] * vector.values[0] + values[1] * vector.values[1] + values[2] * vector.values[2]);	}
	vector4<T> operator* (const T & number) const;
	vector4<T> operator* (const vector4<T> & vector) const;

	private:
		T values[4];
}; 
typedef vector4<float> vec4f;

//--------------------------------------------------------------------------------------------------
// 2D Matrix Structure
//--------------------------------------------------------------------------------------------------

template <typename T>
struct matrix2
{
	matrix2();
	matrix2(vector2<T> x, vector2<T>);

	void operator= (const matrix2 & matrix);

	const T& operator[] (int index) const { return values[index]; }
	T& operator[] (int index) { return values[index]; }

	matrix2<T> operator+ (const matrix2<T> & matrix) const;
	matrix2<T> operator- (const matrix2<T> & matrix) const;
	matrix2<T> operator* (const T & number) const;
	matrix2<T> operator* (const matrix2<T> & input) const;
	vector2<T> operator* (const vector2<T> & input) const;

	matrix2<T> transpose() const;
	T determinant() const;
	matrix2<T> inverse() const;

private:
	T values[4];
};
typedef matrix2<float> mat2f;

//--------------------------------------------------------------------------------------------------
// 3D Matrix Structure
//--------------------------------------------------------------------------------------------------

template <typename T>
struct matrix3
{
	matrix3();
	matrix3(vector3<T> x, vector3<T> y, vector3<T> z);

	void operator= (const matrix3<T> & matrix);

	const T& operator[] (int index) const { return values[index]; }
	T& operator[] (int index) { return values[index]; }

	matrix3<T> operator+ (const matrix3<T> & matrix) const;
	matrix3<T> operator- (const matrix3<T> & matrix) const;
	matrix3<T> operator* (const T & number) const;
	matrix3<T> operator* (const matrix3<T> & input) const;
	vector3<T> operator* (const vector3<T> & input) const;

	matrix3<T> transpose() const;
	T determinant() const;
	matrix3<T> inverse() const;

private:
	T values[9];
};
typedef matrix3<float> mat3f;

//--------------------------------------------------------------------------------------------------
// 4D Matrix Structure
//--------------------------------------------------------------------------------------------------

template <typename T>
struct matrix4
{
	matrix4();
	matrix4(vector4<T> x, vector4<T> y, vector4<T> z, vector4<T> w);

	void operator= (const matrix4<T> & matrix);

	const T& operator[] (int index) const { return values[index]; }
	T& operator[] (int index) { return values[index]; }

	matrix4<T> operator+ (const matrix4<T> & matrix) const;
	matrix4<T> operator- (const matrix4<T> & matrix) const;
	matrix4<T> operator* (const T & number) const;
	matrix4<T> operator* (const matrix4<T> & input) const;
	vector4<T> operator* (const vector4<T> & input) const;

	matrix4<T> transpose() const;
	T determinant() const;
	matrix4<T> inverse() const;

private:
	T values[16];
};
typedef matrix4<float> mat4f;

//--------------------------------------------------------------------------------------------------
// Frustum Structure
//--------------------------------------------------------------------------------------------------

struct frustrum
{
	float fieldOfView;
	float aspectRatio;
	float near;
	float far;
	float top, bottom, left, right;

	frustrum(float FOV = 45.0f, float n = 0.1f, float f = 100.0f, unsigned int imageWidth = 800, unsigned int imageHeight = 600)
	{
		fieldOfView = FOV;
		near = n;
		far = f;
		aspectRatio = (float)imageWidth / (float)imageHeight;

		top = tan(fieldOfView * 0.5 * PI / 180) * near;
		right = top * aspectRatio;
		bottom = top * -1.0f;
		left = right * -1.0f;
	}

	mat4f genProjectionMatrix()
	{
		mat4f temp;

		temp[0] = (2.0f * near) / (right - left);
		temp[1] = 0.0f;
		temp[2] = (right + left) / (right - left);
		temp[3] = 0.0f;

		temp[4] = 0.0f;
		temp[5] = (2.0f * near) / (top - bottom);
		temp[6] = (top + bottom) / (top - bottom);
		temp[7] = 0.0f;

		temp[8] = 0.0f;
		temp[9] = 0.0f;
		temp[10] = -1.0f * (far + near) / (far - near);
		temp[11] = (-2.0 * far * near) / (far - near);

		temp[12] = 0.0f;
		temp[13] = 0.0f;
		temp[14] = -1.0f;
		temp[15] = 0.0f;

		return temp;
	}
};

//--------------------------------------------------------------------------------------------------
// Camera Function:
//
// From	-> To		-> Reference
// Pos	-> Target	-> Up
// Eye	-> Target	-> Up
// Eye	-> Center	-> Up
//--------------------------------------------------------------------------------------------------

mat4f genCameraMatrix(const vec3f & from, const vec3f & to, vec3f reference)
{
	mat4f temp;
	vec3f right, forward, up;

	reference.normalize();

	forward = from - to;
		forward.normalize();

	right = reference * forward;
		right.normalize();

	up = forward * right;

	temp[0] = right[0];		temp[1] = right[1];		temp[2] = right[2];
	temp[4] = up[0];		temp[5] = up[1];		temp[6] = up[2];
	temp[8] = forward[0];	temp[9] = forward[1];	temp[10] = forward[2];

	mat4f trans;
	trans[3] = -from[0];
	trans[7] = -from[1];
	trans[11] = -from[2];

	return (temp * trans);
}

//--------------------------------------------------------------------------------------------------
// Misc
//--------------------------------------------------------------------------------------------------

void print_matrix(mat2f temp)
{
	cout << setprecision(3);
	cout << temp[0] << " " << temp[1] << endl
		<< temp[2] << " " << temp[3] << endl;

	cout << endl;
}

void print_matrix(mat3f temp)
{
	cout << setprecision(3);
	cout << temp[0] << " " << temp[1] << " " << temp[2] << endl
		<< temp[3] << " " << temp[4] << " " << temp[5] << endl
		<< temp[6] << " " << temp[7] << " " << temp[8] << endl;

	cout << endl;
}

void print_matrix(mat4f temp)
{
	cout << setprecision(3);
	cout << temp[0] << " " << temp[1] << " " << temp[2] << " " << temp[3] << endl
		<< temp[4] << " " << temp[5] << " " << temp[6] << " " << temp[7] << endl
		<< temp[8] << " " << temp[9] << " " << temp[10] << " " << temp[11] << endl
		<< temp[12] << " " << temp[13] << " " << temp[14] << " " << temp[15] << endl;

	cout << endl;
}

//--------------------------------------------------------------------------------------------------
// Utilities:
//
// Rotation order from world view:	X, Y, Z
// Rotation order from camera view: -Z, -Y, -X
//
// Transform order: Translation * Rotation * Scale * Vector
//
// MVP: Projection * View * Model
//--------------------------------------------------------------------------------------------------

template <typename T>
matrix4<T> genScalingMatrix(const T & number)
{
	matrix4<T> temp;
		temp[0] = temp[5] = temp[10] = number;
	return temp;
}

template <typename T>
matrix4<T> genScalingMatrix(const T & x, const T & y, const T & z)
{
	matrix4<T> temp;
	temp[0] = x;
	temp[5] = y;
	temp[10] = z;
	return temp;
}

template <typename T>
matrix4<T> genXRotationMatrix(const T & angle)
{
	matrix4<T> temp;
	T rad = angle * (PI / 180.0);

	temp[5] = cos(rad);
	temp[6] = sin(rad);
	temp[9] = sin(rad) * T(-1);
	temp[10] = cos(rad);

	return temp;
}

template <typename T>
matrix4<T> genYRotationMatrix(const T & angle)
{
	matrix4<T> temp;
	T rad = angle * (PI / 180.0);

	temp[0] = cos(rad);
	temp[2] = sin(rad) * T(-1);
	temp[8] = sin(rad);
	temp[10] = cos(rad);

	return temp;
}

template <typename T>
matrix4<T> genZRotationMatrix(const T & angle)
{
	matrix4<T> temp;
	T rad = angle * (PI / 180.0);

	temp[0] = cos(rad);
	temp[1] = sin(rad);
	temp[4] = sin(rad) * T(-1);
	temp[5] = cos(rad);

	return temp;
}


template <typename T>
matrix4<T> genTranslationMatrix(const T & x, const T & y, const T & z)
{
	matrix4<T> temp;
	temp[3] = x;
	temp[7] = y;
	temp[11] = z;
	return temp;
}

//--------------------------------------------------------------------------------------------------
// Main Loop
//--------------------------------------------------------------------------------------------------

int main()
{
	frustrum screen;
	mat4f model, view, projection, mvp;

	print_matrix(model);

	view = genCameraMatrix(vec3f(4, 3, 3), vec3f(0, 0, 0), vec3f(0, 1, 0));

	print_matrix(view);

	projection = screen.genProjectionMatrix();

	print_matrix(projection);

	mvp = projection * view * model;

	print_matrix(mvp);

	std::cin.get();
	std::cin.get();

	return 0;
}

//--------------------------------------------------------------------------------------------------
// 2D Vector Definitions
//--------------------------------------------------------------------------------------------------

template <typename T>
vector2<T>::vector2(T x, T y)
{
	values[0] = x;
	values[1] = y;
}

template <typename T>
void vector2<T>::operator= (const vector2<T> & vector)
{
	values[0] = vector.values[0];
	values[1] = vector.values[1];
}

template <typename T>
vector2<T> vector2<T>::operator+ (const vector2<T> & vector) const
{
	vector2<T> temp;

	temp.values[0] = values[0] + vector.values[0];
	temp.values[1] = values[1] + vector.values[1];

	return temp;
}

template <typename T>
vector2<T> vector2<T>::operator- (const vector2<T> & vector) const
{
	vector2<T> temp;

	temp.values[0] = values[0] - vector.values[0];
	temp.values[1] = values[1] - vector.values[1];

	return temp;
}

template <typename T>
void vector2<T>::normalize()
{
	T len = length();
	if (len > 0)
	{
		len = 1 / len;
		values[0] *= len;
		values[1] *= len;
	}
}

template <typename T>
vector2<T> vector2<T>::operator* (const T & number) const
{
	vector2<T> temp;

	temp.values[0] = values[0] * number;
	temp.values[1] = values[1] * number;

	return temp;
}

template <typename T>
T vector2<T>::operator* (const vector2<T> & vector) const
{
	return(values[0] * vector.values[1] - values[1] * vector.values[0]);
}

//--------------------------------------------------------------------------------------------------
// 3D Vector Definitions
//--------------------------------------------------------------------------------------------------

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
vector3<T> vector3<T>::operator* (const vector3<T> & vector) const
{
	vector3<T> temp;

	temp.values[0] = values[1] * vector.values[2] - values[2] * vector.values[1];
	temp.values[1] = values[2] * vector.values[0] - values[0] * vector.values[2];
	temp.values[2] = values[0] * vector.values[1] - values[1] * vector.values[0];

	return temp;
}

//--------------------------------------------------------------------------------------------------
// 4D Vector Definitions
//--------------------------------------------------------------------------------------------------

template <typename T>
vector4<T>::vector4(T x, T y, T z)
{
	values[0] = x;
	values[1] = y;
	values[2] = z;
	values[3] = 1.0f;
}

template <typename T>
void vector4<T>::operator= (const vector4<T> & vector)
{
	values[0] = vector.values[0];
	values[1] = vector.values[1];
	values[2] = vector.values[2];
}

template <typename T>
vector4<T> vector4<T>::operator+ (const vector4<T> & vector) const
{
	vector4<T> temp;

	temp.values[0] = values[0] + vector.values[0];
	temp.values[1] = values[1] + vector.values[1];
	temp.values[2] = values[2] + vector.values[2];

	return temp;
}

template <typename T>
vector4<T> vector4<T>::operator- (const vector4<T> & vector) const
{
	vector4<T> temp;

	temp.values[0] = values[0] - vector.values[0];
	temp.values[1] = values[1] - vector.values[1];
	temp.values[2] = values[2] - vector.values[2];

	return temp;
}

template <typename T>
void vector4<T>::normalize()
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

template <typename T>
vector4<T> vector4<T>::operator* (const T & number) const
{
	vector4<T> temp;

	temp.values[0] = values[0] * number;
	temp.values[1] = values[1] * number;
	temp.values[2] = values[2] * number;

	return temp;
}

template <typename T>
vector4<T> vector4<T>::operator* (const vector4<T> & vector) const
{
	vector4<T> temp;

	temp.values[0] = values[1] * vector.values[2] - values[2] * vector.values[1];
	temp.values[1] = values[2] * vector.values[0] - values[0] * vector.values[2];
	temp.values[2] = values[0] * vector.values[1] - values[1] * vector.values[0];

	return temp;
}

//--------------------------------------------------------------------------------------------------
// 2D Matrix Definitions
//--------------------------------------------------------------------------------------------------

template <typename T>
matrix2<T>::matrix2()
{
	values[0] = 1;	values[1] = 0;
	values[2] = 0;	values[3] = 1;
}

template <typename T>
matrix2<T>::matrix2(vector2<T> x, vector2<T> y)
{
	values[0] = x[0];	values[1] = x[1];
	values[2] = y[0];	values[3] = y[1];

}

template <typename T>
void matrix2<T>::operator= (const matrix2<T> & matrix)
{
	for (int i = 0; i < 4; i++)
		values[i] = matrix.values[i];
}

template <typename T>
matrix2<T> matrix2<T>::operator+ (const matrix2<T> & matrix) const
{
	matrix2<T> temp;

	for (int i = 0; i < 4; i++)
		temp.values[i] = values[i] + matrix.values[i];

	return temp;
}

template <typename T>
matrix2<T> matrix2<T>::operator- (const matrix2<T> & matrix) const
{
	matrix2<T> temp;

	for (int i = 0; i < 4; i++)
		temp.values[i] = values[i] - matrix.values[i];

	return temp;
}

template <typename T>
matrix2<T> matrix2<T>::operator* (const T & number) const
{
	matrix2<T> temp;

	for (int i = 0; i < 4; i++)
		temp.values[i] = values[i] * number;

	return temp;
}

template <typename T>
matrix2<T> matrix2<T>::operator* (const matrix2<T> & input) const
{
	matrix2<T> temp;

	temp.values[0] = values[0] * input.values[0] + values[1] * input.values[2];
	temp.values[1] = values[0] * input.values[1] + values[1] * input.values[3];
	temp.values[2] = values[2] * input.values[0] + values[3] * input.values[2];
	temp.values[3] = values[2] * input.values[1] + values[3] * input.values[3];

	return temp;
}

template <typename T>
vector2<T> matrix2<T>::operator* (const vector2<T> & input) const
{
	vector2<T> temp;

	temp[0] = values[0] * input[0] + values[1] * input[1];
	temp[1] = values[2] * input[0] + values[3] * input[1];

	return temp;
}

template <typename T>
matrix2<T> matrix2<T>::transpose() const
{
	matrix2<T> temp;

	temp.values[0] = values[0];
	temp.values[1] = values[2];
	temp.values[2] = values[1];
	temp.values[3] = values[3];

	return temp;
}

template <typename T>
T matrix2<T>::determinant() const
{
	return ((values[0] * values[3]) - (values[1] * values[2]));
}

template <typename T>
matrix2<T> matrix2<T>::inverse() const
{
	matrix2<T> temp;

	T det = determinant();
	det = 1 / det;

	temp[0] = values[3];
	temp[1] = values[1] * T(-1);
	temp[2] = values[2] * T(-1);
	temp[3] = values[0];

	temp = temp * det;

	return temp;
}

//--------------------------------------------------------------------------------------------------
// 3D Matrix Definitions
//--------------------------------------------------------------------------------------------------

template <typename T>
matrix3<T>::matrix3()
{
	values[0] = 1;	values[1] = 0;	values[2] = 0;
	values[3] = 0;	values[4] = 1;	values[5] = 0;
	values[6] = 0;	values[7] = 0;	values[8] = 1;

}

template <typename T>
matrix3<T>::matrix3(vector3<T> x, vector3<T> y, vector3<T> z)
{
	values[0] = x[0];	values[1] = x[1];	values[2] = x[2];
	values[3] = y[0];	values[4] = y[1];	values[5] = y[2];
	values[6] = z[0];	values[7] = z[1];	values[8] = z[2];
}

template <typename T>
void matrix3<T>::operator= (const matrix3<T> & matrix)
{
	for (int i = 0; i < 9; i++)
		values[i] = matrix.values[i];
}

template <typename T>
matrix3<T> matrix3<T>::operator+ (const matrix3<T> & matrix) const
{
	matrix3<T> temp;

	for (int i = 0; i < 9; i++)
		temp.values[i] = values[i] + matrix.values[i];

	return temp;
}

template <typename T>
matrix3<T> matrix3<T>::operator- (const matrix3<T> & matrix) const
{
	matrix3<T> temp;

	for (int i = 0; i < 9; i++)
		temp.values[i] = values[i] - matrix.values[i];

	return temp;
}

template <typename T>
matrix3<T> matrix3<T>::operator* (const T & number) const
{
	matrix3<T> temp;

	for (int i = 0; i < 9; i++)
		temp.values[i] = values[i] * number;

	return temp;
}

template <typename T>
matrix3<T> matrix3<T>::operator* (const matrix3<T> & input) const
{
	matrix3<T> temp;

	temp.values[0] = values[0] * input.values[0] + values[1] * input.values[3] + values[2] * input.values[6];
	temp.values[1] = values[0] * input.values[1] + values[1] * input.values[4] + values[2] * input.values[7];
	temp.values[2] = values[0] * input.values[2] + values[1] * input.values[5] + values[2] * input.values[8];

	temp.values[3] = values[3] * input.values[0] + values[4] * input.values[3] + values[5] * input.values[6];
	temp.values[4] = values[3] * input.values[1] + values[4] * input.values[4] + values[5] * input.values[7];
	temp.values[5] = values[3] * input.values[2] + values[4] * input.values[5] + values[5] * input.values[8];

	temp.values[6] = values[6] * input.values[0] + values[7] * input.values[3] + values[8] * input.values[6];
	temp.values[7] = values[6] * input.values[1] + values[7] * input.values[4] + values[8] * input.values[7];
	temp.values[8] = values[6] * input.values[2] + values[7] * input.values[5] + values[8] * input.values[8];

	return temp;
}

template <typename T>
vector3<T> matrix3<T>::operator* (const vector3<T> & input) const
{
	vector3<T> temp;

	temp[0] = values[0] * input[0] + values[1] * input[1] + values[2] * input[2];
	temp[1] = values[3] * input[0] + values[4] * input[1] + values[5] * input[2];
	temp[2] = values[6] * input[0] + values[7] * input[1] + values[8] * input[2];

	return temp;
}

template <typename T>
matrix3<T> matrix3<T>::transpose() const
{
	matrix3<T> temp;

	temp.values[0] = values[0];
	temp.values[1] = values[3];
	temp.values[2] = values[6];
	temp.values[3] = values[1];
	temp.values[4] = values[4]; 
	temp.values[5] = values[7]; 
	temp.values[6] = values[2];
	temp.values[7] = values[5];
	temp.values[8] = values[8];

	return temp;
}

template <typename T>
T matrix3<T>::determinant() const
{
	return ((values[0] * ((values[4] * values[8]) - (values[5] * values[7]))) - (values[1] * ((values[3] * values[8]) - (values[5] * values[6]))) + (values[2] * ((values[3] * values[7]) - (values[4] * values[6]))));
}

template <typename T>
matrix3<T> matrix3<T>::inverse() const
{
	matrix3<T> temp;

	T det = determinant();
	det = 1 / det;

	temp[0] = (values[4] * values[8]) - (values[5] * values[7]);
	temp[1] = ((values[3] * values[8]) - (values[5] * values[6])) * T(-1);
	temp[2] = (values[3] * values[7]) - (values[4] * values[6]);

	temp[3] = ((values[1] * values[8]) - (values[2] * values[7])) * T(-1);
	temp[4] = (values[0] * values[8]) - (values[2] * values[6]);
	temp[5] = ((values[0] * values[7]) - (values[1] * values[6])) * T(-1);

	temp[6] = (values[1] * values[5]) - (values[2] * values[4]);
	temp[7] = ((values[0] * values[5]) - (values[2] * values[3])) * T(-1);
	temp[8] = (values[0] * values[4]) - (values[1] * values[3]);

	temp = temp.transpose();
	temp = temp * det;

	return temp;
}

//--------------------------------------------------------------------------------------------------
// 4D Matrix Definitions
//--------------------------------------------------------------------------------------------------

template <typename T>
matrix4<T>::matrix4()
{
	values[0] = 1;	values[1] = 0;	values[2] = 0;	values[3] = 0;
	values[4] = 0;	values[5] = 1;	values[6] = 0;	values[7] = 0;
	values[8] = 0;	values[9] = 0;	values[10] = 1;	values[11] = 0;
	values[12] = 0;	values[13] = 0;	values[14] = 0;	values[15] = 1;
}

template <typename T>
matrix4<T>::matrix4(vector4<T> x, vector4<T> y, vector4<T> z, vector4<T> w)
{
	values[0] = x[0];	values[1] = x[1];	values[2] = x[2];	values[3] = x[3];
	values[4] = y[0];	values[5] = y[1];	values[6] = y[2];	values[7] = y[3];
	values[8] = z[0];	values[9] = z[1];	values[10] = z[2];	values[11] = z[3];
	values[12] = w[0];	values[13] = w[1];	values[14] = w[2];	values[15] = w[3];
}

template <typename T>
void matrix4<T>::operator= (const matrix4<T> & matrix)
{
	for (int i = 0; i < 16; i++)
		values[i] = matrix.values[i];
}

template <typename T>
matrix4<T> matrix4<T>::operator+ (const matrix4<T> & matrix) const
{
	matrix4<T> temp;

	for (int i = 0; i < 16; i++)
		temp.values[i] = values[i] + matrix.values[i];

	return temp;
}

template <typename T>
matrix4<T> matrix4<T>::operator- (const matrix4<T> & matrix) const
{
	matrix4<T> temp;

	for (int i = 0; i < 16; i++)
		temp.values[i] = values[i] - matrix.values[i];

	return temp;
}

template <typename T>
matrix4<T> matrix4<T>::operator* (const T & number) const
{
	matrix4<T> temp;

	for (int i = 0; i < 16; i++)
		temp.values[i] = values[i] * number;

	return temp;
}

template <typename T>
matrix4<T> matrix4<T>::operator* (const matrix4<T> & input) const
{
	matrix4<T> temp;

	temp.values[0] = values[0] * input.values[0] + values[1] * input.values[4] + values[2] * input.values[8] + values[3] * input.values[12];
	temp.values[1] = values[0] * input.values[1] + values[1] * input.values[5] + values[2] * input.values[9] + values[3] * input.values[13];
	temp.values[2] = values[0] * input.values[2] + values[1] * input.values[6] + values[2] * input.values[10] + values[3] * input.values[14];
	temp.values[3] = values[0] * input.values[3] + values[1] * input.values[7] + values[2] * input.values[11] + values[3] * input.values[15];

	temp.values[4] = values[4] * input.values[0] + values[5] * input.values[4] + values[6] * input.values[8] + values[7] * input.values[12];
	temp.values[5] = values[4] * input.values[1] + values[5] * input.values[5] + values[6] * input.values[9] + values[7] * input.values[13];
	temp.values[6] = values[4] * input.values[2] + values[5] * input.values[6] + values[6] * input.values[10] + values[7] * input.values[14];
	temp.values[7] = values[4] * input.values[3] + values[5] * input.values[7] + values[6] * input.values[11] + values[7] * input.values[15];

	temp.values[8] = values[8] * input.values[0] + values[9] * input.values[4] + values[10] * input.values[8] + values[11] * input.values[12];
	temp.values[9] = values[8] * input.values[1] + values[9] * input.values[5] + values[10] * input.values[9] + values[11] * input.values[13];
	temp.values[10] = values[8] * input.values[2] + values[9] * input.values[6] + values[10] * input.values[10] + values[11] * input.values[14];
	temp.values[11] = values[8] * input.values[3] + values[9] * input.values[7] + values[10] * input.values[11] + values[11] * input.values[15];

	temp.values[12] = values[12] * input.values[0] + values[13] * input.values[4] + values[14] * input.values[8] + values[15] * input.values[12];
	temp.values[13] = values[12] * input.values[1] + values[13] * input.values[5] + values[14] * input.values[9] + values[15] * input.values[13];
	temp.values[14] = values[12] * input.values[2] + values[13] * input.values[6] + values[14] * input.values[10] + values[15] * input.values[14];
	temp.values[15] = values[12] * input.values[3] + values[13] * input.values[7] + values[14] * input.values[11] + values[15] * input.values[15];

	return temp;
}


template <typename T>
vector4<T> matrix4<T>::operator* (const vector4<T> & input) const
{
	vector4<T> temp;

	temp[0] = values[0] * input[0] + values[1] * input[1] + values[2] * input[2] + values[3] * input[3];
	temp[1] = values[4] * input[0] + values[5] * input[1] + values[6] * input[2] + values[7] * input[3];
	temp[2] = values[8] * input[0] + values[9] * input[1] + values[10] * input[2] + values[11] * input[3];
	temp[3] = values[12] * input[0] + values[13] * input[1] + values[14] * input[2] + values[15] * input[3];

	temp[0] /= temp[3];
	temp[1] /= temp[3];
	temp[2] /= temp[3];

	return temp;
}


template <typename T>
matrix4<T> matrix4<T>::transpose() const
{
	matrix4<T> temp;

	temp.values[0] = values[0]; temp.values[1] = values[4]; temp.values[2] = values[8]; temp.values[3] = values[12];
	temp.values[4] = values[1]; temp.values[5] = values[5]; temp.values[6] = values[9]; temp.values[7] = values[13];
	temp.values[8] = values[2]; temp.values[9] = values[6]; temp.values[10] = values[10]; temp.values[11] = values[14];
	temp.values[12] = values[3]; temp.values[13] = values[7]; temp.values[14] = values[11]; temp.values[15] = values[15];

	return temp;
}

template <typename T>
T matrix4<T>::determinant() const
{
	T temp[4];

	temp[0] = (values[5] * ((values[10] * values[15]) - (values[11] * values[14]))) - (values[6] * ((values[9] * values[15]) - (values[11] * values[13]))) + (values[7] * ((values[9] * values[14]) - (values[10] * values[13])));
	temp[1] = (values[4] * ((values[10] * values[15]) - (values[11] * values[14]))) - (values[6] * ((values[8] * values[15]) - (values[11] * values[12]))) + (values[7] * ((values[8] * values[14]) - (values[10] * values[12])));
	temp[2] = (values[4] * ((values[9] * values[15]) - (values[11] * values[13]))) - (values[5] * ((values[8] * values[15]) - (values[11] * values[12]))) + (values[7] * ((values[8] * values[13]) - (values[9] * values[12])));
	temp[3] = (values[4] * ((values[9] * values[14]) - (values[10] * values[13]))) - (values[5] * ((values[8] * values[14]) - (values[10] * values[12]))) + (values[6] * ((values[8] * values[13]) - (values[9] * values[12])));

	return ((values[0] * temp[0]) - (values[1] * temp[1]) + (values[2] * temp[2]) - (values[3] * temp[3]));
}

template <typename T>
matrix4<T> matrix4<T>::inverse() const
{
	matrix4<T> temp;

	T det = determinant();
	det = 1 / det;

	temp[0] = (values[5] * ((values[10] * values[15]) - (values[11] * values[14]))) - (values[6] * ((values[9] * values[15]) - (values[11] * values[13]))) + (values[7] * ((values[9] * values[14]) - (values[10] * values[13])));
	temp[1] = ((values[4] * ((values[10] * values[15]) - (values[11] * values[14]))) - (values[6] * ((values[8] * values[15]) - (values[11] * values[12]))) + (values[7] * ((values[8] * values[14]) - (values[10] * values[12])))) * T(-1);
	temp[2] = (values[4] * ((values[9] * values[15]) - (values[11] * values[13]))) - (values[5] * ((values[8] * values[15]) - (values[11] * values[12]))) + (values[7] * ((values[8] * values[13]) - (values[9] * values[12])));
	temp[3] = ((values[4] * ((values[9] * values[14]) - (values[10] * values[13]))) - (values[5] * ((values[8] * values[14]) - (values[10] * values[12]))) + (values[6] * ((values[8] * values[13]) - (values[9] * values[12])))) * T(-1);

	temp[4] = ((values[1] * ((values[10] * values[15]) - (values[11] * values[14]))) - (values[2] * ((values[9] * values[15]) - (values[11] * values[13]))) + (values[3] * ((values[9] * values[14]) - (values[10] * values[13])))) * T(-1);
	temp[5] = (values[0] * ((values[10] * values[15]) - (values[11] * values[14]))) - (values[2] * ((values[8] * values[15]) - (values[11] * values[12]))) + (values[3] * ((values[8] * values[14]) - (values[10] * values[12])));
	temp[6] = ((values[0] * ((values[9] * values[15]) - (values[11] * values[13]))) - (values[1] * ((values[8] * values[15]) - (values[11] * values[12]))) + (values[3] * ((values[8] * values[13]) - (values[9] * values[12])))) * T(-1);
	temp[7] = (values[0] * ((values[9] * values[14]) - (values[10] * values[13]))) - (values[1] * ((values[8] * values[14]) - (values[10] * values[12]))) + (values[2] * ((values[8] * values[13]) - (values[9] * values[12])));

	temp[8] = (values[1] * ((values[6] * values[15]) - (values[7] * values[14]))) - (values[2] * ((values[5] * values[15]) - (values[7] * values[13]))) + (values[3] * ((values[5] * values[14]) - (values[6] * values[13])));
	temp[9] = ((values[0] * ((values[6] * values[15]) - (values[7] * values[14]))) - (values[2] * ((values[4] * values[15]) - (values[7] * values[12]))) + (values[3] * ((values[4] * values[14]) - (values[6] * values[12])))) * T(-1);
	temp[10] = (values[0] * ((values[5] * values[15]) - (values[7] * values[13]))) - (values[1] * ((values[4] * values[15]) - (values[7] * values[12]))) + (values[3] * ((values[4] * values[13]) - (values[5] * values[12])));
	temp[11] = ((values[0] * ((values[5] * values[14]) - (values[6] * values[13]))) - (values[1] * ((values[4] * values[14]) - (values[6] * values[12]))) + (values[2] * ((values[4] * values[13]) - (values[5] * values[12])))) * T(-1);

	temp[12] = ((values[1] * ((values[6] * values[11]) - (values[7] * values[10]))) - (values[2] * ((values[5] * values[11]) - (values[7] * values[9]))) + (values[3] * ((values[5] * values[10]) - (values[6] * values[9])))) * T(-1);
	temp[13] = (values[0] * ((values[6] * values[11]) - (values[7] * values[10]))) - (values[2] * ((values[4] * values[11]) - (values[7] * values[8]))) + (values[3] * ((values[4] * values[10]) - (values[6] * values[8])));
	temp[14] = ((values[0] * ((values[5] * values[11]) - (values[7] * values[9]))) - (values[1] * ((values[4] * values[11]) - (values[7] * values[8]))) + (values[3] * ((values[4] * values[9]) - (values[5] * values[8])))) * T(-1);
	temp[15] = (values[0] * ((values[5] * values[10]) - (values[6] * values[9]))) - (values[1] * ((values[4] * values[10]) - (values[6] * values[8]))) + (values[2] * ((values[4] * values[9]) - (values[5] * values[8])));

	temp = temp.transpose();
	temp = temp * det;

	return temp;
}