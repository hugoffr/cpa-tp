#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>

using namespace std;

#define SYSTEMTIME clock_t

void f0(bool *primes, long long k, long long n)
{
    do
    {
        for (long long j = k*k ; j<n ; j+=2*k)
        {   
            primes[j>>1]=true;
        }
        
        // K goes up to the next, skipping one because even numbers are never prime
        do
        {
            k+=2;
        } while (k*k <= n && primes[k>>1]);
        
    } while (k*k <= n);
}

void f1(bool *primes, long long k, long long n)
{

}

void f2(bool *primes, long long k, long long n)
{

}

void f3(bool *primes, long long k, long long n)
{
    
}

int main (int argc, char *argv[])
{

    SYSTEMTIME Time1, Time2;

    long long n;
    
    cout << "Power of 2: ";
    cin >> n;

    int op = 0;

    do
    {
        cout << endl << "Which function method be used?" << endl;
        cout << "0. Premade (Teacher's Code)" << endl;
        cout << "1. First" << endl;
        cout << "2. Second" << endl;
        cout << "3. Third" << endl;
        cin >> op;

    } while (op < 0 && op > 3); // Don't forget to update these if you add options
 
    n = pow(2,n);
    
    bool *primes = new bool[n];
    
    long long k = 3;

    Time1 = clock();
    switch (op)
    {
        case 0:
            f0(primes, k, n);
            break;
        case 1:
            f1(primes, k, n);
            break;
        case 2:
            f2(primes, k, n);
            break;
        case 3:
            f3(primes, k, n);
            break;
        default:
            cout << endl << "Something went wrong." << endl;
            break;
    }

    Time2 = clock();

    char st[100];
    sprintf(st, "Time: %3.3f seconds\n", (double)(Time2 - Time1) / CLOCKS_PER_SEC);
    cout << st;
    
    /*
    for (int i=1; i<n; i+=2)
        if (!primes[i>>1])
            cout << i << " ";
    */
}