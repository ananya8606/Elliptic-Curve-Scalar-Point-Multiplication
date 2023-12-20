#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>

typedef struct ECC_Point {
    int X, Y;
} ECC_Point;

typedef struct EllipticCurve {
    int A, B, p; // Coefficients and prime field
} EllipticCurve;

// int modInverse(int a, int m) {
//     int m0 = m;
//     int y = 0, x = 1;

//     if (m == 1) {
//         return 0;
//     }

//     while (a > 1) {
//         int q = a / m;
//         int t = m;

//         m = a % m, a = t;
//         t = y;

//         y = x - q * y;
//         x = t;
//     }

//     if (x < 0) {
//         x += m0;
//     }

//     return x;
// }

int power(int base, int exp, int mod) {
    int result = 1;
    base = base % mod;
    while (exp > 0) {
        if (exp % 2 == 1)
            result = (result * base) % mod;
        exp = exp >> 1;
        base = (base * base) % mod;
    }
    return result;
}

// modular inverse using Fermat Little a^(p-1) = 1modp
int modInverse(int a, int p) {
    return power(a, p - 2, p);
}

ECC_Point addPoints(ECC_Point P, ECC_Point Q, EllipticCurve curve) {
    if (P.X == Q.X && (P.Y + Q.Y) % curve.p == 0) {
        ECC_Point infinity = {INT_MIN, INT_MIN};
        return infinity;
    }

    int m;
    if (P.X != Q.X || P.Y != Q.Y) {
        m = ((Q.Y - P.Y + curve.p) * modInverse(Q.X - P.X + curve.p, curve.p)) % curve.p;
    } else {
        m = ((3 * P.X * P.X + curve.A) * modInverse(2 * P.Y, curve.p)) % curve.p;
    }

    int X3 = (m * m - P.X - Q.X) % curve.p;
    if (X3 < 0)
        X3 = curve.p + X3;
    int Y3 = (m * (P.X - X3) - P.Y) % curve.p;
    if (Y3 < 0)
        Y3 = curve.p + Y3;
    ECC_Point result = {X3, Y3};
    //printf("add: %d %d\n", result.X, result.Y);
    return result;
}

ECC_Point doublePoint(ECC_Point P, EllipticCurve curve) {
    if (P.X == INT_MIN && P.Y == INT_MIN) {
        return P;
    }

    int m = ((3 * P.X * P.X + curve.A) * modInverse(2 * P.Y, curve.p)) % curve.p;
    int X3 = (m * m - 2 * P.X) % curve.p;
    if (X3 < 0)
        X3 = curve.p + X3;
    int Y3 = (m * (P.X - X3) - P.Y) % curve.p;
    if (Y3 < 0)
        Y3 = curve.p + Y3;
    ECC_Point result = {X3, Y3};
    //printf("doubling: %d %d\n", result.X, result.Y);
    return result;
}

ECC_Point scalarMultiplication(int k, ECC_Point P, EllipticCurve curve) {
    ECC_Point Q = {INT_MIN, INT_MIN};
    ECC_Point G = P;

    for (int i = 31; i >= 0; i--) {
        Q = doublePoint(Q, curve);
        if (((k >> i) & 1) == 1) {
            if (Q.X == INT_MIN && Q.Y == INT_MIN) {
                Q = G;
            } else {
                Q = addPoints(Q, G, curve);
            }
        }
    }

    return Q;
}



int isPointOnCurve(ECC_Point P, EllipticCurve curve) {
    int left = (P.Y * P.Y) % curve.p;
    int right = (P.X * P.X * P.X + curve.A * P.X + curve.B) % curve.p;
    return (left == right);
}

int main() {
    EllipticCurve curve;
   // int a, b, p;
    printf("Enter valid (a, b, p): ");
    scanf("%d %d %d", &curve.A, &curve.B, &curve.p);
    //EllipticCurve curve = {a,b,p};

    ECC_Point basePoint;
    printf("Enter P coordinates (x y): ");
    scanf("%d %d", &basePoint.X, &basePoint.Y);

    int n;
    printf("Enter the scalar value for point multiplication: ");
    scanf("%d", &n);



    //int privateKey = 21;
    if (!isPointOnCurve(basePoint, curve)) {
        printf("Base point is not on the elliptic curve.\n");
        return 1;
    }
    clock_t start = clock();
    ECC_Point result = scalarMultiplication(n, basePoint, curve);
    clock_t end = clock();

    double time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Final Result: (%d, %d)\n", result.X, result.Y);
    printf("Execution time: %f seconds\n", time);
    return 0;
}
