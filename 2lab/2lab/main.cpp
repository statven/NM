#include "newton_method.h"
#include <iostream>

using namespace std;

int main() {
    cout << "M = 0.01:" << endl;
    newton_method(1, 1, df1_dx1, df1_dx2, df2_dx1, df2_dx2);
    cout << endl;

    cout << "M = 0.05:" << endl;
    newton_method(1, 1, df1_dx1, df1_dx2, df2_dx1, df2_dx2);
    cout << endl;

    cout << "M = 0.1:" << endl;
    newton_method(1, 1, df1_dx1, df1_dx2, df2_dx1, df2_dx2);

    return 0;
}