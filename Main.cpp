#include <iostream>
#include <string>
#include "VectorMatrix.h"
#include "Algo.h"
using namespace std;

void readData(string filename, vvi& A, vi& y0) {
    ifstream fin(filename);

    int n, m;
    fin >> n >> m;
    A.resize(n);
    for (auto& a : A) a.resize(m);
    fin >> A;

    fin >> n;
    y0.resize(n);
    fin >> y0;

    fin.close();
}

void getReport(string filename, vd error) {
    ofstream fout(filename);
    fout << error;
    fout.close();
}

void doTest(string testName, int wmin = 0) {
    vvi A;
    vi y0;
    vd error;
    vi x;

    readData(testName + ".txt", A, y0);
    cout << testName << endl;

    for (int j = 1; j <= 3; j++) {
        cout << "Algo" << j << endl;
        switch (j) {
        case 1: {
            int iter = Algo1(A, y0, x, error, wmin);
            getReport(testName + "_algo" + to_string(j) + ".txt", error);
            cout << "Iteration number:" << iter << endl;
            cout << x;
            break; }
        case 2: {
            int iter = Algo2(A, y0, x, error, wmin);
            getReport(testName + "_algo" + to_string(j) + ".txt", error);
            cout << "Iteration number:" << iter << endl;
            cout << x;
            break; }
        case 3: {
            int iter = Algo3(A, y0, x, error, wmin);
            getReport(testName + "_algo" + to_string(j) + ".txt", error);
            cout << "Iteration number:" << iter << endl;
            cout << x;
            break; }
        default:
            break;
        }

    }
    cout << endl;
}

int main() {
    doTest("test1", 5);
    doTest("test2_1", 5);
    doTest("test2_2", 5);
    doTest("test3", 5);
    doTest("test4_1", 5);
    doTest("test4_2", 5);
    doTest("test4_3", 5);
    doTest("test4_4", 5);
}