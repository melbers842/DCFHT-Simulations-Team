//
// Created by Angel Madrigal on 2/5/26.
//
#include <iostream>
#include "MeshTest.h"

using std::cout;
using std::endl;

int main() {
    Eigen::MatrixXd m(2, 2);
    m(0, 0) = 3;
    m(0, 1) = 2.5;
    m(1, 0) = -1;
    m(1, 1) = m(0, 1) + m(1,0);


    cout << m << endl;

    return 0;
}
