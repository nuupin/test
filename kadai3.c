#include <stdio.h>
#include <math.h>

int main(void) {
    double f, f1, f2, dx, x;

    f1 = 1;
    dx = 0.01;

    for (int i=0; i<=100; ++i) {
        x = i*dx;
        f2 = f1 + dx*x*f1;

        f = exp(x*x/2);

        printf("x = %f 数値解 -> f1 = %f 正確な値 -> f = %f \n", x, f1, f);

        f1 = f2;
    }
}
