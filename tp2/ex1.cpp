#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <omp.h>

using namespace std;

#define SYSTEMTIME clock_t

void f0(bool *primes, long long n)
{
    long long k = 3;

    do
    {
        for (long long j = k*k ; j<n ; j+=2*k)
        {   
            primes[j>>1]=true;
        }
        
        do
        {
            k+=2;
        } while (k*k <= n && primes[k>>1]);
        
    } while (k*k <= n);
}

void f1(bool *primes, long long n)
{
    long long k = 3;

    do
    {
        // Mark all multiples of k between k^2 and n
        // by dividing each element by k and checking if the ramainder is 0
        for (long long j = k*k; j<=n; j++)
        {
            if (j % k == 0) 
            {
                primes[j] = true;
            }
        }
        
        // Smallest unmarked number becomes K
        do
        {
            k+=2;
        } while (k*k <= n && primes[k]);
        
    } while (k*k <= n); 
}

void f2(bool *primes, long long n)
{
    long long k = 3;
    do
    {
        // Mark all multiples of k between k^2 and n
        // by fast-marking all values that are computed
        for (long long j = k*k; j<=n; j+=k*2)
        {
                primes[j>>1] = true;
        }
        
        // Smallest unmarked number becomes K
        do
        {
            k+=2;
        } while (k*k <= n && primes[k]);
        
    } while (k*k <= n); 
}

void f3(bool *primes, long long n)
{
    long long k = 3;
    int cacheBlockSize = 65536; //my L1 cache has 256kiB, we are setting this to 64kiB
    long long seedList[(int)(sqrt(n)+1)] = {}; //List to store the seeds
    long long foundSeeds = 0; //How many seeds have been found so far

    long long blockNumber = n / cacheBlockSize + 1; //In how many blocks we will divide data
    bool allSeedsFound = false; //Wether all seeds have been found or not

    for(long long block = 0; block < blockNumber; block++)
    {
        //First test seeds already found
        long long blockStart = block*cacheBlockSize;
        long long blockEnd = (block+1)*cacheBlockSize;
        if (blockEnd > n)
        {
            blockEnd = n+1;
        }

        for(long long seed_i = 0; seed_i < foundSeeds; seed_i++)
        {
            long long seed = seedList[seed_i];
            //Lets say we are in block starting in 100 000 and we are testing the seed 7. We can start at 0 and add 7 until we reach 100 000 or we can divide 100 000 by 7,
            //round up and multiply by 7 and we have the first multiple of 7 within the current block.

            //Second implementation
            long long markingStart = blockStart;
            long long squaredSeed = seed*seed;
            if (blockStart < squaredSeed) //We must start testing from seed*seed, so if the block start before we must start a little after the beggining
            {
                if(squaredSeed >= blockEnd)
                {
                    break; //If the square of the seed is not inside this block than no other square of seed is going to be inside, so we move on
                }
                markingStart = squaredSeed;
            }

            long long j = (int)(ceil(blockStart/((double)seed)) * seed);

            while (j < blockEnd)
            {
                primes[j] = true;
                j+=seed;
            }
        }

        if(allSeedsFound) //Check if there are still seeds to find
        {
            continue;
        }

        long long k = blockStart;

        if(k == 0)
        {
            k = 3;
        }

        while(primes[k] && k < blockEnd) //Find first unmarked number
        {
            k++;
        }

        if (k == blockEnd) //Case in which all numbers are already marked in this block
        {
            continue;
        }

        do
        {
            seedList[foundSeeds] = k;
            foundSeeds++;
            // Mark all multiples of k between k^2 and n
            // by fast-marking all values that are computed
            long long j = k*k;
            while(j < blockEnd)
            {
                primes[j] = true;
                j+=k;
            }
            
            // Smallest unmarked number becomes K
            do
            {
                k+=2;
            } while (k*k <= n && primes[k]);
            
        } while (k*k < blockEnd);

        if (k*k <= n)
        {
            allSeedsFound = true;
        } 
    }
}

void f4(bool *primes, long long n)
{
    long long k = 3;
    int cacheBlockSize = 65536; //my L1 cache has 256kiB, we are setting this to 64kiB
    long long seedList[(int)(sqrt(n)+1)] = {}; //List to store the seeds
    long long foundSeeds = 0; //How many seeds have been found so far

    long long blockNumber = n / cacheBlockSize + 1; //In how many blocks we will divide data
    bool allSeedsFound = false; //Wether all seeds have been found or not

    #pragma omp parallel for
    for(long long block = 0; block < blockNumber; block++)
    {
        //First test seeds already found
        long long blockStart = block*cacheBlockSize;
        long long blockEnd = (block+1)*cacheBlockSize;
        if (blockEnd > n)
        {
            blockEnd = n+1;
        }

        for(long long seed_i = 0; seed_i < foundSeeds; seed_i++)
        {
            long long seed = seedList[seed_i];
            //Lets say we are in block starting in 100 000 and we are testing the seed 7. We can start at 0 and add 7 until we reach 100 000 or we can divide 100 000 by 7,
            //round up and multiply by 7 and we have the first multiple of 7 within the current block.

            //Second implementation
            long long markingStart = blockStart;
            long long squaredSeed = seed*seed;
            if (blockStart < squaredSeed) //We must start testing from seed*seed, so if the block start before we must start a little after the beggining
            {
                if(squaredSeed >= blockEnd)
                {
                    break; //If the square of the seed is not inside this block than no other square of seed is going to be inside, so we move on
                }
                markingStart = squaredSeed;
            }

            long long j = (int)(ceil(blockStart/((double)seed)) * seed);

            while (j < blockEnd)
            {
                primes[j] = true;
                j+=seed;
            }
        }

        if(allSeedsFound) //Check if there are still seeds to find
        {
            continue;
        }

        long long k = blockStart;

        if(k == 0)
        {
            k = 3;
        }

        while(primes[k] && k < blockEnd) //Find first unmarked number
        {
            k++;
        }

        if (k == blockEnd) //Case in which all numbers are already marked in this block
        {
            continue;
        }

        do
        {
            seedList[foundSeeds] = k;
            foundSeeds++;
            // Mark all multiples of k between k^2 and n
            // by fast-marking all values that are computed
            long long j = k*k;
            while(j < blockEnd)
            {
                primes[j] = true;
                j+=k;
            }
            
            // Smallest unmarked number becomes K
            do
            {
                k+=2;
            } while (k*k <= n && primes[k]);
            
        } while (k*k < blockEnd);
        
        if (k*k <= n)
        {
            allSeedsFound = true;
        } 
    }
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
        cout << endl << "Which function method should be used?" << endl;
        cout << "0. Premade (Teacher's Code)" << endl;
        cout << "1. Normal Method" << endl;
        cout << "2. Fast Marking" << endl;
        cout << "3. Segmented Fast Marking" << endl;
        cout << "4. Multi-Core Segmented Fast Marking" << endl;
        cin >> op;

    } while (op < 0 && op > 3); // Don't forget to update these if you add options
 
    n = pow(2,n);
    
    bool *primes = new bool[n];

    Time1 = clock();

    switch (op)
    {
        case 0:
            f0(primes, n);
            break;

        case 1:
            f1(primes, n);
            break;

        case 2:
            f2(primes, n);
            break;

        case 3:
            f3(primes, n);
            break;

        case 4:
            f4(primes, n);
            break;

        default:
            cout << endl << "Something went wrong." << endl;
            break;
    }

    Time2 = clock();

    char st[100];
    sprintf(st, "Time: %3.3f seconds\n", (double)(Time2 - Time1) / CLOCKS_PER_SEC);
    cout << st;
    
    if(false){
        if (op == 0)
        {
            for (int i=1; i<n; i+=2)
                if (!primes[i>>1])
                    cout << i << " ";
        }
        else
        {
            for (int i=1; i<n; i+=2)
                if (!primes[i])
                    cout << i << " ";
        }
    }
}