#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <iomanip>

using namespace std;

vector<vector<double>> matrix_mult(vector<vector<double>> A, vector<vector<double>> B) {
    vector<vector<double>> C(A.size(), vector<double>(B[0].size(), 0));
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < B[i].size(); ++j) {
            for (int z = 0; z < A.size(); ++z) {
                C[i][j] += A[i][z] * B[z][j];
            }
        }
    }

    return C;
}

vector<double> matrix_mult(vector<vector<double>> A, vector<double> B) {
    vector<double> C(B.size(), 0);
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A.size(); ++j) {
            C[i] += A[i][j] * B[j];
        }
    }

    return C;
}

vector<vector<double>> make_p_matrix(int n, int swap1 = 0, int swap2 = 0) {
    vector<vector<double>> res(n, vector<double>(n, 0));
    for (int i = 0; i < n; ++i) { res[i][i] = 1; }
    swap(res[swap1], res[swap2]);

    return res;
}


vector<vector<double>> make_m_matrix(vector<vector<double>> matr, int iter) {
    vector<vector<double>> res(matr.size(), vector<double>(matr.size(), 0));
    for (int i = 0; i < res.size(); ++i) { res[i][i] = 1; }
    res[iter][iter] = 1 / matr[iter][iter];

    for (int i = iter + 1; i < res.size(); ++i) {
        res[i][iter] = -(matr[i][iter] / matr[iter][iter]);
    }

    return res;
}


vector<double> matrix_gauss_method(vector<vector<double>> matr, vector<double> ans) {
    vector<vector<double>> P;
    vector<vector<double>> M;

    for (int j = 0; j < matr.size(); ++j) {
        int max = 0, n = 0;
        for (int i = j; i < matr.size(); ++i) {
            if (fabs(matr[i][j]) > max) {
                max = fabs(matr[i][j]);
                n = i;
            }
        }

        P = make_p_matrix(matr.size(), j, n);

        matr = matrix_mult(P, matr);
        ans = matrix_mult(P, ans);



        M = make_m_matrix(matr, j);

        matr = matrix_mult(M, matr);
        ans = matrix_mult(M, ans);

    }


    for (int i = matr.size() - 2; i >= 0; i--) {
        for (int j = i + 1; j < matr.size(); j++) {
            ans[i] -= matr[i][j] * ans[j];
        }
    }

    return ans;
}



void print_matrix(vector<vector<double>> matr) {
    cout << '\n';
    for (int i = 0; i < matr.size(); i++) {
        cout << "\t|";
        for (int j = 0; j < matr[i].size(); j++) {
            cout << setw(7) << setprecision(3) << matr[i][j] << "|";
        }
        cout << '\n';
    }
}



double dichotomy_method(long double (*func)(double[], double[], int, double), double a[], double b[], int k, double left, double right, double eps) {
    if (func(a, b, k, left) * func(a, b, k, right) > 0.0) {
        throw std::logic_error("Wrong arguments");
    }

    eps = fabs(eps);

    while (right - left > eps) {
        double middle = (left + right) / 2.0;
        if (func(a, b, k, middle) * func(a, b, k, right) < 0.0) {
            left = middle;
        }
        else {
            right = middle;
        }
    }
    return (left + right) / 2.0;
}

double func(double x) {
    return pow(M_E, -x) + x * x - 2;
}



double* cheb_zeros(double a, double b, int n) {
    double* zeros = new double[n];
    cout << "Chebyshev polynom zeros\n\t";
    for (int i = 0; i < n; ++i) {
        zeros[i] = ((a + b) / 2) + ((b - a) / 2) * cos(((2 * i + 1) * M_PI) / (2 * n));
        cout << zeros[i] << ' ';
    }   cout << '\n';
    return zeros;
}

double* generate(double (*func)(double), double*& a, double*& b, int n, double from, double to) {
    a = new double[n];
    b = new double[n];
    double* zeros = cheb_zeros(from, to, n);

    for (int i = 0; i < n; ++i) {
        a[i] = zeros[n - i - 1];
        b[i] = func(a[i]);
        // from += zeros[i];
    }

    double* h = new double[n];
    for (int i = 1; i < n; ++i) {
        h[i] = a[i] - a[i - 1];
    }
    return h;
}


long double lag(double a[], double b[], int k, double x) {
    long double ans = 0;

    for (int i = 0; i < k; ++i) {
        long double a1 = 1;

        for (int j = 0; j < k; ++j) {
            if (j == i) continue;
            a1 *= (x - a[j]) / (a[i] - a[j]);
        }

        ans += a1 * b[i] ;
    }

    return ans;
}




void cubic_spline(double* h, int n, double x[], double y[]) {

    cout << fixed << setprecision(3) << "\n\n\n\tSpline building\n\n";
    vector<vector<double>> A(n - 2, vector<double>(n - 2, 0));
    for (int i = 0; i < n - 2; ++i) {
        for (int j = 0; j < n - 2; ++j) {
            if (i == j) {
                A[i][j] = (h[i + 2] + h[i + 1]) / 3;
            }

            if (i == j - 1) {
                A[i][j] = h[j + 1] / 6;
            }

            if (i - 1 == j) {
                A[i][j] = h[i + 1] / 6;
            }
        }
    }

    vector<vector<double>> H(n - 2, vector<double>(n, 0));

    for (int i = 0; i < n - 2; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                H[i][j] = 1 / h[i + 1];
            }

            if (i == j - 1) {
                H[i][j] = (-1 / h[i + 1]) - (1 / h[j + 1]);
            }

            if (i == j - 2) {
                H[i][j] = 1 / h[i + 2];
            }

        }
    }
    cout << "Matrix A: \n";
    print_matrix(A);
    cout << "Matrix H: \n";
    print_matrix(H);

    vector<double> f(n);
    for (int i = 0; i < n; ++i) f[i] = y[i];

    vector<double> Bf = matrix_mult(H, f);

    vector<double> c = matrix_gauss_method(A, Bf);

    for (int i = n - 2; i > 0; --i) c[i] = c[i - 1];
    c[0] = 0;
    c[n - 1] = 0;
    cout << "Vector c: \n";
    for (auto i : c) cout << i << ' ';
    cout << '\n';

    cout << "Splines: \n";
    for (int i = 1; i < n; ++i) {
        double a = y[i];
        double d = (c[i] - c[i - 1]) / h[i];
        double b = ((h[i] / 2) * c[i]) - ((h[i] * h[i] * d) / 6) + ((y[i] - y[i - 1]) / h[i]);
        cout << i << ": " << a << " + " << b << "(x - " << x[i] << ") + " <<
            c[i] / 2 << "(x - " << x[i] << ")^2 + " <<
            d / 6 << "(x - " << x[i] << ")^3" <<
            "   x is [" << x[i - 1] << "; " << x[i] << "]\n";

    }

}

int main()
{
    double* x;
    double* y;
    double from = -1;
    double to = 0;
    int n = 10;
    double* h = generate(func, x, y, n, from, to);
    cout << "x and y:\n";
    for (int i = 0; i < n; ++i) {
        cout << "i: " << i << "\tx: " << x[i] << "\ty: " << y[i] << '\n';
    }
    cout << "h (step): ";
    for (int i = 1; i < n; ++i) {
        cout << h[i] << ", ";
    }   cout << '\n';

    cout << "Method 1, using interpolation:\n" << dichotomy_method(lag, x, y, n, from, to, 0.0001) << '\n';
    cout << "Method 2, using reversed interpolation\n" << lag(y, x, n, 0.0) << '\n';

    cubic_spline(h, n, x, y);
    int asd = 0;
    cin >> asd;


    return 0;
}