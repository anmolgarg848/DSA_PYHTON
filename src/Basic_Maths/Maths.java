package Basic_Maths;

import java.util.List;

import static java.lang.Math.log;
import static java.lang.Math.min;

public class Maths {
    public static void main(String[] args) {
        dig_extr(222333666);
        System.out.println();
        System.out.println(rev_num(25678945));
        System.out.println(is_pal(232));
        System.out.println(is_armstrng(1634));
        prnt_All_div1(36);
        System.out.println();
        System.out.println(check_prime(20));
        System.out.println(GCD_HCF(20, 40));
        System.out.println(euclidean_algo(20, 40));
    }

    //    EXTRACTION OF DIGITS
//    T.C --> log10(n) if the number of iterations is based on division
//    the time complexity will always be logarithmic
    public static void dig_extr(int n) {
        while (n > 0) {
            int dig = n % 10;
            System.out.print(dig + " ");
            n = n / 10;
        }
    }

    //    T.C - O(1) constant
    public static int cnt_dig(int n) {
        if (n == 0) return 1; // Special case for 0
        int cnt = (int) Math.log10(Math.abs(n)) + 1;
        return cnt;
    }

    public static int rev_num(int n) {
        int rn = 0;
        while (n > 0) {
            int dig = n % 10;
            rn = (rn * 10) + dig;
            n = n / 10;
        }
        return rn;
    }

    public static boolean is_pal(int n) {
        return n == rev_num(n);
    }

    public static boolean is_armstrng(int n) {
        int sum = 0;
        int x = n;
        int tot_dig = cnt_dig(n);
        while (x > 0) {
            int dig = x % 10;
            sum += Math.pow(dig, tot_dig);
            x = x / 10;
        }
        return n == sum;
    }

    //    print all divisors
//   TC-> O(n)
    public static void prnt_All_div(int n) {
        for (int i = 1; i <= n; i++) {
            if (n % i == 0) {
                System.out.print(i + " ");
            }
        }

    }

    //T.C -> O sqrt(n)
    public static void prnt_All_div1(int n) {
        //sqrt fn is called again and takes some time so uses i*i
        //        for (int i = 1; i<=Math.sqrt(n); i++) {
        for (int i = 1; i * i <= n; i++) {
            if (n % i == 0) {
                System.out.print(i + ", ");
                if ((n / i) != i) {
                    System.out.print(n / i + ", ");
                }
            }
        }

    }

    //Prime number is a number which exactly has 2 factors (1 & itself)
//T.C- Osqrt(n)
    public static boolean check_prime(int n) {
        if (n <= 1) {
            return false; // 0 and 1 are not prime numbers
        }
        int cnt = 0;
        for (int i = 1; i * i <= n; i++) {
            if (n % i == 0) {  //here 6 is already count as a factor
                cnt++;
//example count  2 and 18 both as divisor of 36 at same time
//                condition is applied not to count 6,6 two times ex-> 36/6=6
                if (n / i != i) { // If n/i is not equal to i, count both divisors
                    cnt++;
                }
            }
        }
        return cnt == 2; // A prime number has exactly 2 divisors
    }

    /*PRINT ALL PRIME FACTORS OF NUMBER
    BRUTE FORCE --> loop till 1 to n if its divisor then check for prime or loop til root of n
    T.C--> O(ROOT N  * (2 * ROOT N) )
OPTIMIZED --> factorize with prime only start from 2
    */

    public static void print_prime_factors(int n) {
        //4 can't be divisor by default as we have already extended 2 so any multiple of 2 can't be
        //divisor
        for (int i = 2; i <= n; i++) {
            if (n % i == 0) {
                System.out.println(i + " ");
                while (n % i == 0) {
                    n = n / i;
                }
            }
        }
        //T.C  --> O(N) example 37 for a large number that is prime

    }

    public static void print_prime_factors_optimized(int n) {
        //4 can't be divisor by default as we have already extended 2 so any multiple of 2 can't be divisor
        //
        for (int i = 2; i <= (int) Math.sqrt(n); i++) {
            if (n % i == 0) {
                System.out.println(i + " ");
                while (n % i == 0) { //n is changing here so i*i will not work
                    n = n / i;
                }
            }
        }
        if (n != 1) System.out.println(n); //only prime number left after divisor's
        //T.C  --> O(root N* log n) example 37 for a large number that is prime

    }

    //PRINT PRIME NUMBERS TILL N --> SIEVE OF ERADTHONESES
    public static void sieve(int n) {
        //T.C --> O(LOG(LOG N) ) + O(N)
        int[] box = new int[n + 1];
        for (int i = 2; i*i<= n; i++) {
            if (box[i] == 0) {
                    //multiple of prime cannot be prime as well
                    for (int j = i * i; j<=n ; j+=i) { //shift of i
                        box[j]= 1;
                    }
            }
        }
        for (int i = 2; i <= n; i++) {
            if (box[i] == 0) {
                System.out.println(i + " ");
            }
        }

    }

    //POWER EXPONENTIATION
    public static int power_expo(int x, int n) {
        int ans = 1;
        //odd even approach
        while (n < 0) {
            if (n % 2 == 0) {
                n = n / 2;
                x = x * x;
            } else {//odd
                n = n - 1;
                ans *= x;
            }
        }
        //if n is -ve
        if (n < 0) return 1 / ans;
//        if n is double just make the ans double and divide by 1.0 if its -ve

        return ans;
    }


    //      GCD /HCF greatest common divisor/highest common factor
//T.C -> O(min(x,y))
    public static int GCD_HCF(int x, int y) {
        int n = 1;
        for (int i = 1; i <= Math.min(x, y); i++) {
            if (x % i == 0 && y % i == 0) {
                n = i;
            }
        }
        /*
//     apprch 2--> Reverse loop
//     works for lot of cases like (20,40) but in case of (11,13) still
//        whole loop will run till 1
        for (int i = min(x,y); i >=1 ; i--) {
            if (x%i==0&&y%i==0){
                n = i;
                break;
            }
        }

         */
        return n;
    }

//    Euclidean Algorithm

    public static int euclidean_algo(int a, int b) {
/*ALGORITHM
gcd(a,b) = gcd(a-b,b)        a>b always
gcd(20,15) = gcd(5,15)
gcd(15,5) = gcd(10,5)
gcd(10,5) = gcd(5,5)
gcd(5,5) = gcd(0,5)   till  a or b one of them gets zero

5 is the answer
Observation above  subtraction gives 5 times iteration leds to higher tc
we can achieve same using  modulo in doubt look for ex -> gcd(52,10)
gcd(a,b) = gcd(a%b,b)        a>b always
gcd(20,15) = gcd(20%15,15)
gcd(5,15) = gcd(15%5,5)
gcd(0,5) = gcd(5,0)
5 is the answer in less iteration
T.C-> O(log(min(a,b))

 */
        while (a > 0 && b > 0) {
            if (a > b) {
                a = a % b;
            } else {
                b = b % a;
            }
        }
        if (a == 0) {
            return b;
        }
        return a;
    }
}