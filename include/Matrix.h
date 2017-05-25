#ifndef MATRIX_CLASS
#define MATRIX_CLASS

#include "Complex.h"
#include <cstdlib>
#include <cmath>
//#include "State.h"

const double L_PI = M_PI;

class Matrix {
	private:
		int rows;
		int cols;
		Complex ***data;
		bool notamatrix;
		//void reshape(int,int);
	public:
		//Constructors and destructors
		Matrix(int inRows, int inCols);
		Matrix(const Matrix&);
		//Matrix(const State*);
		~Matrix();

		//Members functions
		Complex E(int,int) const;						//Get value at given indices
		Complex& E(int,int);
		void E(int,int,Complex);						//Assign value at given indices with Complex number
		void E(int,int,double);							//Assign value at given indices with single double. Real part remains unchanged
		void E(int,int,double,double);					//Assign value at given indices with two doubles. Both real and imaginary parts are assigned.
		void E(int,int,bool,double);					//Assign value at given indices with effectively one double. Assigns only imaginary part and bool is ignored. (By convention use true)
		int getRows() const;
		int getCols() const;
		bool isSquare();
		int isVec() const;								//Returns -1 for row vec, 0 for not a vec, 1 for col vec
		double magnitude() const;							//Returns length if matrix is a vector (row or col)
		int length() const;								//Returns number of elements on main axis if vector
		bool multCompatible(Matrix&);
		Matrix multiply(Matrix&);
		//Vector multiply(Vector&);
		Matrix inverse();
		Matrix transpose();
		Matrix adjoint();
		Complex determinant();
		bool isNAM() const;
		Matrix add(Matrix& inM);
		Matrix tensorProduct(const Matrix&);
		static Matrix nTensorM(const Matrix&, int inN);
		void clone(const Matrix& inM);

		//Display functions
		void print() const;

		//Operator overloads
			//Assignment
		void operator=(const Matrix&);
			//Addition and subtraction
		Matrix operator+(Matrix&);
		Matrix operator-(Matrix&);
			//Scalar and matrix multiplication
		Matrix operator*(double);						//Note: This means scalar multiplication always has to be at the END of a statement involving matrices
		Matrix operator*(Matrix&);
		//Vector operator*(Vector&);
			//Tensor product
		Matrix operator^(const Matrix&);
		Matrix operator^=(const Matrix&);
			//Vector element access
		Complex operator[](int index) const;
		Complex& operator[](int index);

		//Matrix generators
		static Matrix colVec(int);
		static Matrix rowVec(int);
		static Matrix identity(int);
		static Matrix identity();
		static Matrix X();
		static Matrix Y();
		static Matrix Z();
		static Matrix H();
		static Matrix R(double theta);
		static Matrix CNOT();
		static Matrix CMAT(const Matrix& inMat);
		static Matrix X(int Nbit,int nth);		//All nths or bit indices start from 1 in this context
		static Matrix Y(int Nbit,int nth);
		static Matrix Z(int Nbit,int nth);
		static Matrix H(int Nbit,int nth);
		static Matrix R(int Nbit,int nth, double theta);
		static Matrix O(int Nbit,int nth);
		static Matrix CNOT(int Nbit, int cbit, int obit);
		static Matrix CMAT(const Matrix& inMat, int Nbit, int cbit, int obit);
		static Matrix PERM(int Nbits, int cbit, int Mbits, int inA, int inC);
        static Matrix CNOT(int Nbit, int* qargs, int nargs = -1);//qargs[0] = target, quargs[1 -> (n-1)] = controls
        static Matrix CS(int Nbit, int c1, int c2);
        static Matrix CCZ(int Nbit, int c1, int c2, int c3);

		static Matrix loadMatrix(const char* filename);

		//Shorthand
		Matrix T();
};

#endif
