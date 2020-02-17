#include <iostream>
#include <bits/stdc++.h>
#include <math.h>
#include <Eigen/Dense>
#include "Simplex.h"

using namespace std;
using namespace Eigen;

int main(){

    int m = 3;
    int n = 6;
    int a = 1;
    int u = 2;
    
    MatrixXd mat(m, n);

    mat << -5, -2, 0, 0, 0, 0,
            10, 2, 1, 0, 0, 60,
            10, 0, 0, -1, 1, 30;

    Simplex solver(mat, m, n, a, u);
    if(solver.hasSolution){
	    vector<pair<int, double>> solution = solver.solution;
	    vector<pair<int, double>>::iterator it;
        for(it = solution.begin(); it < solution.end(); it++){
		    cout << "X" << it->first << " = " << it->second << endl;
	    }
        cout << "Optimal : " << solver.optimal << endl;
    }
	return 0;
}
