#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

class Matrix
{
private:
    int m;
    int n;
public:
    double** matrix;
    Matrix(int m, int n)
    {
        this->m = m;
        this->n = n;
        matrix = new double* [m];
        for (int i = 0; i < m; i++)
            matrix[i] = new double[n];
    }
    void set()
    {
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                cin >> matrix[i][j];
    }
    int getM()
    {
        return m;
    };
    int getN()
    {
        return n;
    };
    void print()
    {
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (j == n - 1) cout << setprecision(4) << fixed<< matrix[i][j];
                else cout << setprecision(4) << fixed << matrix[i][j] << " ";
            }
            cout << endl;
        }
    }
    friend double** operator + (Matrix& a, Matrix& b)
    {
        int answerM = a.getM();
        int answerN = b.getN();
        Matrix answer(answerM, answerN);
        for (int i = 0; i < answerM; i++)
        {
            for (int j = 0; j < answerN; j++)
            {
                answer.matrix[i][j] = a.matrix[i][j] + b.matrix[i][j];
            }
        }
        return answer.matrix;
    }
    friend double** operator - (Matrix& a, Matrix& b)
    {
        int answerM = a.getM();
        int answerN = b.getN();
        Matrix answer(answerM, answerN);
        for (int i = 0; i < answerM; i++)
        {
            for (int j = 0; j < answerN; j++)
            {
                answer.matrix[i][j] = a.matrix[i][j] - b.matrix[i][j];
            }
        }
        return answer.matrix;
    }
    friend double** operator * (Matrix& a, Matrix& b)
    {
        int answerM = a.getM();
        int answerN = b.getN();
        int help = a.getN();
        Matrix answer(answerM, answerN);
        for(int i = 0; i < answerM; i++)
            for(int j = 0; j < answerN; j++) {
                answer.matrix[i][j] = 0;
                for (int k = 0; k < help; ++k) {
                    answer.matrix[i][j] += (a.matrix[i][k] * b.matrix[k][j]);
                }
            }
        return answer.matrix;
    }
    double** transpose(Matrix& a)
    {
        int answerN = a.getN();
        int answerM = a.getM();
        Matrix answer(answerN, answerM);
        for (int i = 0; i < answerN; i++)
        {
            for (int j = 0; j < answerM; j++)
            {
                answer.matrix[i][j] = a.matrix[j][i];
            }
        }
        return answer.matrix;
    }
    double** operator = (Matrix& b)
    {
        Matrix temp (b.m, b.n);
        temp.matrix = b.matrix;
        return temp.matrix;
    }
};
class SquareMatrix : public Matrix
{
public:
    SquareMatrix(int n) : Matrix(n, n) {}
    SquareMatrix(double** matrix, int n) : Matrix(n, n) {}
};
class EliminationMatrix : public SquareMatrix {
public:
    EliminationMatrix(Matrix a, int n, int i, int j): SquareMatrix(n) {
        for (int i = 0; i < n; i++) {
            matrix[i][i] = 1;
        }
        matrix[i][j] = -(a.matrix[i][j]/a.matrix[j][j]);
    }
};
class PermutationMatrix : public SquareMatrix {
public:
    PermutationMatrix(int n, int i, int j): SquareMatrix(n) {
        for (int k = 0; k < n; k++) {
            if (k == i) {
                matrix[k][j] = 1;
            } else if (k == j) {
                matrix[k][i] = 1;
            } else {
                matrix[k][k] = 1;
            }
        }
    }
};
class IdentityMatrix : public Matrix {
public:
    IdentityMatrix(int n): Matrix(n, n) {
        for (int i = 0; i < n; i++) {
            matrix[i][i] = 1;
        }
    }
    double** inverse(Matrix& a)
    {
        int answerN = a.getN();
        IdentityMatrix answer(answerN);
        for (int i = 0; i < answerN; i++)
            for (int i = 0; i < answerN - 1; i++) {
                //permutations
                double maxi = abs(a.matrix[i][i]);
                int indJ = i;
                for (int j = i+1; j < answerN; j++) {
                    if (maxi < abs(a.matrix[j][i])) {
                        maxi = abs(a.matrix[j][i]);
                        indJ = j;
                    }
                }
                if (indJ != i) {
                    PermutationMatrix p(answerN, indJ, i);
                    a.matrix = p * a;
                    answer.matrix = p * answer;
                }
                //elimination to get U matrix
                for (int j = i + 1; j < answerN; j++) {
                    if (a.matrix[j][i] != 0) {
                        EliminationMatrix e(a, answerN, j, i);
                        a.matrix = e * a;
                        answer.matrix = e * answer;
                    }
                }
            }
        //elimination to get D matrix
        for (int i = answerN - 1; i > 0; i--) {
            for (int j = i - 1; j >= 0; j--) {
                if (a.matrix[j][i] != 0) {
                    EliminationMatrix e(a, answerN, j, i);
                    a.matrix = e * a;
                    answer.matrix = e * answer;
                }
            }
        }
        //Diagonal normalization
        for (int i = 0; i < answerN; i++) {
            for (int j = 0; j < answerN; j++) {
                answer.matrix[i][j] /= a.matrix[i][i];

            }
        }
        return answer.matrix;
    }
};
class ColumnVector : public Matrix {
public:
    ColumnVector(int n): Matrix(n, 1) {}
    double& operator[](int i) {
        return matrix[i][0];
    }
};
bool equal_or_not(Matrix& a, Matrix& b, double eps)
{
    if (a.getM() != b.getM() || a.getN() != b.getN())
    {
        return false;
    }

    for (int i = 0; i < a.getN(); i++)
    {
        for (int j = 0; j < a.getM(); j++)
        {
            if (abs(a.matrix[i][j] - b.matrix[i][j]) > eps)
            {
                return false;
            }
        }
    }

    return true; // Matrices are equal.
}
ostream& operator <<(ostream& out, Matrix& m)
{
    m.print();
    return (out);
}

istream& operator >>(istream& in, Matrix& m)
{
    m.set();
    return (in);
}

void InteriorPointAlgorithm(Matrix A, Matrix C, Matrix x_sol_new, double eps, double alpha, int nA, int mA, int nC)
{
    ColumnVector x_sol_old(nA);

    SquareMatrix P(nA);

    IdentityMatrix D(nA); // diagonal matrix of trial solution


    Matrix A_new(mA, nA);
    ColumnVector C_new(nC);

    Matrix A_new_transpose(nA, mA);
    SquareMatrix AAt(mA);
    IdentityMatrix AAt_inv(mA);
    Matrix AtAAti(nA, mA);
    Matrix AtAAtiA(nA, nA);
    IdentityMatrix I(nA);

    ColumnVector ones(nA);
    for(int i = 0; i < nA; i++) {
        ones.matrix[i][0] = 1;
    }
    
    bool flag;
    ColumnVector cp(nC);
    while(!(equal_or_not(x_sol_new, x_sol_old, eps))) {
        for(int i = 0; i < nA; i++) {
            D.matrix[i][i] = x_sol_new.matrix[i][0];
        }
        x_sol_old.matrix = x_sol_new.matrix;
        //cout<<A;
        //cout<<C;
        A_new.matrix = A * D;

        C_new.matrix = D * C;

        // calculations for finding P
        A_new_transpose.matrix = A_new.transpose(A_new);
        AAt.matrix = A_new * A_new_transpose;
        AAt_inv.matrix = AAt.matrix;
        AAt_inv.matrix = AAt_inv.inverse(AAt_inv);
        AtAAti.matrix = A_new_transpose * AAt_inv;
        AtAAtiA.matrix = AtAAti * A_new;
        P.matrix = I - AtAAtiA;
        cp.matrix = P * C_new;
        // we need to find min neg elem from cp min_neg
        double min_neg = cp.matrix[0][0];
        flag=false;
        if (min_neg < 0){
            flag = true;
        }

        for(int i = 1; i < nC; i++) {
            if (cp.matrix[i][0] < 0){
                flag = true;
            }
            if (cp.matrix[i][0] < min_neg && cp.matrix[i][0] < 0) {
                min_neg = cp.matrix[i][0];
            }
        }
        if(flag==false){
            break;
        }
        
        double v = abs(min_neg);
        ColumnVector x_temp(nA);
        ColumnVector avc_temp(nA);
        for(int i = 0; i < nA; i++) {
            avc_temp.matrix[i][0] = cp.matrix[i][0] * alpha / v;
        }
        x_temp.matrix = ones + avc_temp;
        x_sol_new.matrix = D * x_temp;
    }
    if(flag==false){
        cout<<"The problem does not have solution!\n";
    }
    else{
        cout << x_sol_new;
    }
}

int main() {
    int nC;
    cout << "Enter the number of the rows of the C: \n" ;
    cin >> nC;
    ColumnVector C(nC); // A vector of coefficients of objective function - C.
    cout << "Enter the C: \n" ;
    cin >> C;

    int mA, nA;
    cout << "Enter the number of the rows of the A: \n" ;
    cout << "Enter the number of the columns of the A: \n" ;
    cin >> mA >> nA;
    Matrix A(mA, nA); // A matrix of coefficients of constraint function - A.
    cout << "Enter the A: \n" ;
    cin >> A;

    int nB;
    cin >> nB;
    cout << "Enter the number of the rows of  the b: \n" ;
    ColumnVector b(nB); // A vector of right-hand side numbers - b.
    cout << "Enter the b: \n" ;
    cin >> b;

    double eps;
    cout << "Enter the the epsilon: \n" ;
    cin >> eps; // The approximation accuracy ϵ.
    double alpha1 = 0.5;
    double alpha2 = 0.9;
    ColumnVector x_sol_new(nA); // trial solution x
    cout<< "Enter the initial trial solution: \n";
    cin >> x_sol_new;

    cout << "For alpha = 0.5:\n";
    InteriorPointAlgorithm(A,C, x_sol_new, eps, alpha1, nA, mA, nC);
    cout << "For alpha = 0.9:\n";
    InteriorPointAlgorithm(A,C, x_sol_new, eps, alpha2, nA, mA, nC);
}