package Arrays;

import java.util.*;

public class Arrays_Medium {
    public static void main(String[] args) {
        int[] arr = {1, 2, 3, 1, 1, 1, 1, 2, 5};
//        int[] arr = {1, 2, 3};
        int k = 10;
//        Gen_sub_Brute1(arr, k);
//        Gen_sub_Brute2(arr, k);
        pascal_optimal(6);
    }

    //Brute force -> Generate all subarrays
    //Longest subarray with sum equals to k
    public static void Gen_sub_Brute1(int[] arr, int a) {
//        T.C --> O(n^3)
        int n = arr.length;
        int len = 0;
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                int sum = 0;
                for (int k = i; k <= j; k++) {
//                    System.out.print(arr[k] + " ");
                    sum += arr[k];
                }
                if (sum == a) {
                    len = Math.max(len, (j - i + 1));
                }
//                System.out.println();
            }
        }
        System.out.println(len);
    }

    //    Brute 2 --> Reducing unnecessary loop based on observation
    public static void Gen_sub_Brute2(int[] arr, int a) {
//        T.C --> O(n^2)
        int n = arr.length;
        int len = 0;
        for (int i = 0; i < n; i++) {
            int sum = 0;
            for (int j = i; j < n; j++) {
                sum += arr[j];
                if (sum == a) {
                    len = Math.max(len, (j - i + 1));
                }
            }
        }
        System.out.println(len);
    }

    public static int Gen_sub_Better(int[] arr, int k) {

//    Approach -->  Hashing + PrefixSum
//        it also handles if the array has +ves and -ves
//        Logic> to get k there must be x-k should be there then only rest subarray is sums up to k
//       T.C--> O(N LOG N)     S.C --> O(N)
//Worst case if lot of collisions happens    T.C--> O(N^N)
        HashMap<Integer, Integer> map = new HashMap<>(); //stores sum , index
        int preSum = 0, max_len = 0;
        map.put(0, -1);
        for (int i = 0; i < arr.length; i++) {
            int ele = arr[i];
            preSum += ele;
            if (map.containsKey(preSum - k)) {
                max_len = Math.max(max_len, i - map.get(preSum - k));
            }
            // Ensure we only store the first occurrence to maximize the subarray length
            if (!map.containsKey(preSum)) {
                map.put(preSum, i);
            }
        }
        return max_len;
    }

    //Optimal --> only if the array contains +ves and zeroes
    public static int Gen_sub_optimal(int[] arr, int k) {
//    Algo --> Greedy /two pointer
//    T.C --> O(2N) worst case  S.C-> O(1)
        int start = 0, end = 0;
        int sum = 0, maxLen = 0;
        while (end < arr.length) {
            sum += arr[end];
            // Shrink the window as long as the sum is greater than k
            while (sum > k && start <= end) { //this loop is not running for n every time for sure
//Since both pointers independently traverse the array at most n times,
// the total operations are O(n + n), which simplifies to O(n).
                sum -= arr[start];
                start++;
            }
            // If the current sum equals k, update the max length
            if (sum == k) {
                maxLen = Math.max(maxLen, end - start + 1);
            }
            end++;
        }
        return maxLen;
    }

    //cnt all subarrays whose sum equals to K
//    Approach --> If you have X looking for k then you have to remove Those elements that sums up to
//    X-K
    public static int cnt_subarrays_k(int[] arr, int k) {
//STORES (prefix_sum, count)
        //T.C--> O(N * LOG N)
        //S.C --> O(N)
        int n = arr.length;
        Map<Integer, Integer> map = new HashMap<>();
        map.put(0, 1);
        int preSum = 0;
        int cnt = 0;
        for (int i = 0; i < n; i++) {
            preSum += arr[i];
            if (map.containsKey(preSum - k)) {  //worst case log N if collision happens
                cnt += map.get(preSum - k);
            }
            map.put(preSum, map.getOrDefault(preSum, 0) + 1);

        }

        return cnt;
    }

    //TWO SUM --> Two varieties
    public int[] twoSum(int[] arr, int tar) {
//   T.C --> O(N)  -- when map works for o(1)  best and average case of unoredred map in worse
//    case when collision happens O(n^2)
//   T.C --> O(N LOG N)  -- when map works for o(logn)
//S.C--> O(n)  using map dumping every element


        Map<Integer, Integer> mp = new HashMap<>();
        int[] ans = {-1, -1};
        for (int i = 0; i < arr.length; i++) {
            int ele = arr[i];
            if (mp.containsKey(tar - ele)) {
                return new int[]{mp.get(tar - ele), i};
            } else {
                mp.put(ele, i);
            }
        }

        return ans;
    }

    public int[] twoSum_sort(int[] arr, int tar) {
//    Approach 2--> Binary Search
//saves space only for yes or no get not optimal for variety 2  as indices to be preserved need
//        sorting  + datastructure then traversal
//     (ele,idx)   [(1,2),(3,0),(8,1)]

        return new int[]{-1, -1};
    }


    /*Sort 0's, 1's and 2's
    Brute -- use merge sort
    Better --   count and overwrite array T.C--> o(2N)   using no extra space just modifying the given array
    Optimal->
    */
    public static void sort_zeroes(int[] arr) {
//Dutch national flag algorithm/2 pointers

//        000000000000 ....         1111111111     Unsorted .......  2222222 rightmost part
//       [0 --> low-1]leftmost part [low ->mid-1] [mid ---> high]    [high+1 --> n-1]


        int lo = 0, n = arr.length;
        int mid = 0;
        int hi = n - 1;
        while (mid <= hi) {
            if (arr[mid] == 0) {
                swap(arr, mid, lo);
                lo++;
                mid++; //shrink
            } else if (arr[mid] == 1) {
                mid++;
            } else { //mid ==2
                swap(arr, mid, hi);
                hi--;
            }

        }
//         time complexity of O(n) and space complexity of O(1)
    }

    public static void swap(int[] arr, int i, int j) {
        int tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }
//    Majority element --> ele appears more than N/2
    /*
    BRUTE --> TWO LOOPS   o(n^2)
    BETTER --> HASHING    O( NLOGN ) + O(N) ,  S.C -->O(N)  ORDERED MAP USED
    OPTIMAL --> MOORE'S VOTING ALGORITHM  INTUITION --> ELEMENT'S CANCEL'S OUT EACH OTHER
    LOGIC --> for element appear's more than n/2 times , other element requires same cnt to cancel it out
     */

    public static int maj_ele(int[] arr) {
        int n = arr.length;
        int ele = arr[0], cnt = 1;
//        find element
        for (int i = 1; i < n; i++) {
            if (cnt == 0) {
                ele = arr[i];
                cnt++;
            } else {
                if (ele == arr[i]) cnt++;
                else cnt--;
            }

        }
//        verify--> THIS STEP NOT REQUIRED WHEN PROBLEM STATES ARRAY ALWAYS HAVE MAJ_ELEMENT
        cnt = 0;
        for (int i = 0; i < n; i++) {
            if (arr[i] == ele) {
                cnt++;
            }
        }
//T.C --> O(N)
        return cnt > (n / 2) ? ele : -1;
    }


    //Kadane's Algorithm --> Maximum Subarray sum
    public static int Kadanes_Algo(int[] arr) {
//Intuition --> consider only those elements who contributes to get maximum sum of subarray
        int pre_sum = 0, n = arr.length, max = Integer.MIN_VALUE;
        for (int i = 0; i < n; i++) {
            int ele = arr[i];
            if (pre_sum < 0) {  //reset as -ves will not contribute to get maximum sum
                pre_sum = 0;
            }
            pre_sum += ele;
            max = Math.max(max, pre_sum);
        }
//to consider empty subarray
//        if (max<0)max = 0;

        return max;
    }

    //Print subarray which gives maximum sum
//    T.C--> O(n) S.C --> O(1)
    public static void print_subarray(int[] arr) {
//   maintain j and i every time max is updated as soon as sum reset's means new subarray starts
//Intuition --> consider only those elements who contributes to get maximum sum of subarray
        int pre_sum = 0, n = arr.length, max = Integer.MIN_VALUE;
        int ansStarts = -1, ansEnds = -1, potential_start = -1;
        for (int i = 0; i < n; i++) {
            int ele = arr[i];
            if (pre_sum == 0) {
                potential_start = i; //new potential subarray starts
            }

            if (pre_sum < 0) {  //reset as -ves will not contribute to get maximum sum
                pre_sum = 0;
            }
            pre_sum += ele;
            if (max < pre_sum) {
                max = pre_sum;
                ansStarts = potential_start;
                ansEnds = i;
            }


        }
//to consider empty subarray
//        if (max<0)max = 0;
        for (int i = ansStarts; i <= ansEnds; i++) {
            System.out.println(arr[i] + " ");
        }
    }

    //    Best time to buy & sell stocks
    public static int BuySellStocks(int[] arr) {
//        buying at a low price and selling at a higher price later.
        int buy = arr[0], prft = 0;
        for (int ele : arr) {
            if (ele < buy) {   //buy it
                buy = ele;
            } else { //sell it
                prft = Math.max(prft, ele - buy);
            }

        }
        return prft;
    }

    //    Rearrange array elements by sign

    //Brute Force --> divide array into +ve and -ve and merge
    public int[] rearrangeArray(int[] arr) {
//        Intuition --> +ves need to place at even and -ves at odd
//        T.C --> O(N + N/2) , S.C --> O(N/2)
        int n = arr.length;
        int[] pos = new int[n / 2];
        int[] neg = new int[n / 2];
        int p = 0, m = 0;
        for (int i = 0; i < n; i++) {
            int ele = arr[i];
            if (ele > 0) {
                pos[p++] = ele;
            } else {
                neg[m++] = ele;
            }
        }
//fill elements
        for (int i = 0; i < (n / 2); i++) {
            arr[2 * i] = pos[i]; //even
            arr[2 * i + 1] = neg[i];//odd
        }


        return arr;
    }

    //    Reducing two pass into one pass
    public int[] rearrangeArray_Better(int[] arr) {
        int n = arr.length;
        int[] ans = new int[n];
        int pos = 0, neg = 0;
        for (int i = 0; i < n; i++) {
            int ele = arr[i];
            if (ele < 0) { //-ve
                ans[2 * neg + 1] = ele;
                neg++;
            } else { //+ve
                ans[2 * pos] = ele;
                pos++;
            }
        }
        return ans;
    }
//    Optimal Inplace --> try later

    public int[] rearrangeArray_Optimal(int[] arr) {
        int n = arr.length;
        /* ALGO
 Find the first positive at an odd index
 Find the first negative  at an even index until both is <n and swap
//        Try Later
        */


        return arr;
    }

    //    Rearrange Variety 2 --> numbers ara not equally divided if any +ves and -ves left add them
//    in the end without altering the order
//    Approach --> Fall back to brute force solution of variety 1;

    public int[] rearrangeArray_Variety2(int[] arr) {
//Case --> pos!=-ve   can be pos>-ve ,or  == or , pos<-ve
//        Try Later in revision Session
//        T.C --> O(2N) , SC -> O(N)
        return arr;
    }


    /*
Total permutations --> n!
Different Ways to Generate Permutations of an Array:
->Simple Recursive Algorithm for Generating All Permutations of an Array
->Iterative Algorithm for Generating All Permutations of an Array
->Heap’s Algorithm for Generating All Permutations of an Array(Recommended)
->Steinhaus Johnson Trotter Algorithm for Generating All Permutations of an Array


 */
    //    Genarate all Permutations
    public static void Generate_permutations(int[] arr, List<Integer> ds, List<List<Integer>> ans, boolean[] freq) {
 /*
T.C --> N! (GENERATING N FACTORIAL ANS) * N ( FOR LOOPING TILL N )
S.C--> O(N) FOR DATA STRUCTURE  + O(N) FOR VISITED ARRAY/HASH ARRAY
Auxillary space (depth of recursion ) --> O(N)
  */
        if (ds.size() == arr.length) {
            ans.add(new ArrayList<>(ds));
            return;
        }
        for (int i = 0; i < arr.length; i++) {
            if (!freq[i]) {
                freq[i] = true;
                Generate_permutations(arr, ds, ans, freq);
                ds.remove(ds.size() - 1);
                freq[i] = false;
            }
        }

    }

    //    Driver method for all permutation approach
    public static List<List<Integer>> call_permute(int[] arr) {
        List<List<Integer>> ans = new ArrayList<>();
        boolean[] freq = new boolean[arr.length];
        List<Integer> ds = new ArrayList<>();
        Generate_permutations(arr, ds, ans, freq); //Brute force
        Gen_perm_Optimized(arr, 0, ans);//optimized
        return ans;
    }

    //    Optimized Approach --> T.C --> O(N!) * O(N)  S.C --> O(1)   , A.S -->O(N)
    public static void Gen_perm_Optimized(int[] arr, int idx, List<List<Integer>> ans) {
        if (idx == arr.length) {
            List<Integer> tmp = new ArrayList<>();
            for (int ele : arr) {
                tmp.add(ele);
            }
            ans.add(new ArrayList<>(tmp));
            return;
        }
        for (int i = idx; i < arr.length; i++) {
            swap(arr, idx, i);
            Gen_perm_Optimized(arr, idx + 1, ans);
            swap(arr, idx, i);
        }
    }


    //     PERMUTATION
    public static void Next_Permutation_Brute(int[] arr) {
/*
1 . Genarate all permutation in sorted Order.
2. Linear Search for current arr.
3. return next_index  if not available then fall back to start using modulo.

   NOTE-> just tell this approach thoroughly as its waste of time to implement it as it takes way much higher T.C --> N!*N
   example --> for array of length 15 it will take almost 10^12 iterations which is impossible
    */
        /*INTERVIEW TIP ALWAYS FOLLOWS

        OBSERVATION --> ALGORITHM --> DRY_RUN --> CODE

         */

    }

    public static void Next_Permutation_Optimal(int[] arr) {
    /*INTUITION
//    [2 1 5 4 3 0 0]
    1. find longest prefix match   (THE BREAKPOINT)    ( raj < rax < rbx )
    2. find some ele >than the break point element but the next possible closest/smallest
    3.Try to place remaining in sorted order
     */
//T.C --> O(3 N )
        int n = arr.length;
        int idx = -1;
//        step 1--> The break point
        for (int i = n - 2; i >= 0; i--) {
            if (arr[i] < arr[i + 1]) {
                idx = i;
                break;  //found the longest prefix
            }
        }
//        edge case
        if (idx == -1) {
            rev(arr, 0, n - 1);
            return;
        }

//    step 2 --> get the next smallest
        for (int i = n - 1; i > idx; i--) {
            if (arr[i] > arr[idx]) {
                swap(arr, i, idx);
                break;
            }
        }
//step 3--> Try to place remaining in sorted order to get the next smallest
        rev(arr, idx + 1, n - 1);

    }

    public static void rev(int[] arr, int i, int j) {
        while (i < j) {
            swap(arr, i, j);
            i++;
            j--;
        }
    }

    //    Leaders in an array
//    Brute  --> Loop on each element O(N^2)
    public static List<Integer> leaders_optimal(int[] arr) {
        List<Integer> ans = new ArrayList<>();
        int cur_max = -1;
        for (int i = arr.length - 1; i >= 0; i--) {
            int ele = arr[i];
            if (ele > cur_max) {
                cur_max = ele;
                ans.add(ele);
            }
        }

        return ans;
    }

    //Longest Consecutive Sequence
    public static int Long_cons_seq_Brute(int[] arr) {
//       T.C --> O(N^3)
        int len = 1, n = arr.length;
        for (int i = 0; i < n; i++) {
            int cnt = 1;
            int ele = arr[i];
            while (linearSearch(arr, ele) == true) {
                cnt++;
                ele++;
            }
            len = Math.max(cnt, len);
        }
        return len;
    }

    //    Linear search
    public static boolean linearSearch(int[] arr, int target) {
        int n = arr.length;
        for (int i = 0; i < n; i++) {
            if (arr[i] == target) {
                return true;
            }
        }
        return false;
    }

    //    Sort the array
    public static int Long_cons_seq_Better(int[] arr) {
//        T.C--> O(NLOG N) + O(N)
        Arrays.sort(arr);
        int n = arr.length, len = 1, cnt = 0, last_smaller = Integer.MIN_VALUE;
        for (int i = 0; i < n; i++) {
            if (arr[i] - 1 == last_smaller) {
                cnt++;
                last_smaller = arr[i];
            } else if (arr[i] != last_smaller) { //reset and skip duplicates
                cnt = 1;
                last_smaller = arr[i];
            }
            len = Math.max(len, cnt);

        }
        return len;
    }


    //    Using SET
    public static int Long_cons_seq_Optimal(int[] arr) {
        if (arr.length < 1) return 0;
        Set<Integer> st = new HashSet<>();
        for (int ele : arr) {   //O(N)
            st.add(ele);  //O(1) ASSUME
        }
        int len = 1;
        for (int ele : st) {
            if (!st.contains(ele - 1)) {//if previous is not available means its a starting point
                int cnt = 0;
                int cur = ele;
                while (st.contains(cur)) { //hardly take 2N iteration at max
                    cnt++;
                    cur = cur + 1;
                }
                len = Math.max(len, cnt);
            }

        }
        return len;
//        T.C --> N + 2N  = 3N Assuming  set insertion is O(1) not logN  otherwise brute is better as it increases space complexity
//        S.C --> O(N)

    }

    //    SET MATRIX ZEROES
    /*Brute  --> Two pass
    1st pass --> mark -1
    2nd pass --> mark -1 as zeroes apart from zeroes
        T.C--> (N*M) * (N+M) + (n*m)
     */
    public static void setZeroes_brute(int[][] matrix) {
        int n = matrix.length, m = matrix[0].length;
        for (int r = 0; r < n; r++) {
            for (int c = 0; c < m; c++) {
                if (matrix[r][c] == 0) {
                    //mark row
                    for (int i = 0; i < n; i++) {
                        if (matrix[r][i] != 0) {
                            matrix[r][i] = -1;
                        }
                    }
                    //mark col
                    for (int i = 0; i < m; i++) {
                        if (matrix[i][c] != 0) {
                            matrix[i][c] = -1;
                        }
                    }
                }
            }
        }
//2nd pass
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                if (matrix[i][j] == -1) {
                    matrix[i][j] = 0;
                }
            }
        }
    }

    public static void setZeroes_better(int[][] matrix) {
//        T.C -> 2* N*M , S.C--> O(N) + O(M)

        int n = matrix.length, m = matrix[0].length;
        int[] ro = new int[n];
        int[] col = new int[m];

        for (int r = 0; r < n; r++) {
            for (int c = 0; c < m; c++) {
                if (matrix[r][c] == 0) {
                    ro[r] = 1;
                    col[c] = 1;
                }
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                if (ro[i] == 1 || col[j] == 1) {
                    matrix[i][j] = 0;
                }
            }
        }

    }

    //    Approach keep a track in the matrix itself
    public static void setZeroes_optimal(int[][] matrix) {
//        T.C--> 2* (N*M)
        int n = matrix.length, m = matrix[0].length;
        // col    matrix[0][...]
        // row    matrix[..][0]
        int col0 = 1;
        //mark
        for (int r = 0; r < n; r++) {
            for (int c = 0; c < m; c++) {
                if (matrix[r][c] == 0) {
                    matrix[r][0] = 0;
                    if (c != 0) {
                        matrix[0][c] = 0;
                    } else {
                        col0 = 0;
                    }
                }
            }
        }

        for (int i = 1; i < n; i++) {
            for (int j = 1; j < m; j++) {
                if (matrix[i][0] == 0 || matrix[0][j] == 0) {
                    matrix[i][j] = 0;
                }
            }
        }
        //mark 0th row and 0th col
        if (matrix[0][0] == 0) {
            for (int i = 0; i < m; i++) {
                matrix[0][i] = 0;
            }
        }
        if (col0 == 0) {
            for (int i = 0; i < n; i++) {
                matrix[i][0] = 0;
            }
        }
    }

    //ROTATE MATRIX IMAGE BY 90 DEGREES
    //Brute Force --> put elements at its correct position
    public static void Rotate_image(int[][] arr) {
        int m = arr[0].length, n = arr.length;
        int[][] ans = new int[n][m];
//        decoded that --> [i][j] => [j][n-1-i]

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                ans[j][n - 1 - i] = arr[i][j];
            }
        }
    }

    //T.C --> O(N^2) , S.C =>O(N^2)
    //OPTIMAL --> TRANSPOSE + REVERSE
    /*HOW TO TRANSPOSE
    diagonals will stay remain same (0,0),(1,1),(2,2) ....
     */
    public static void Rotate_image_optimal(int[][] arr) {

//        T.C --> O(N/2 * N/2)
        transpose(arr);

//        O(N * N/2(TWO POINTER APPROACH TO REVERSE)
        for (int[] a : arr) {

            rev(a); //RUNS FOR N/2 TIMES
        }
        //now the array is rotated

    }

    public static void transpose(int[][] arr) {
//transpose means rows will become columns and vice-versa
        int m = arr[0].length, n = arr.length;
        //traverse only right half
        /*
        i      i+1  j
        0 --> {1 to n-1}
        1 --> {2 to n-1}

        i---> 0 to n-2
        j--> i+1 to n-1
         */


        for (int i = 0; i <= n - 2; i++) {
            for (int j = i + 1; j <= n - 1; j++) {
                int tmp = arr[i][j];
                arr[i][j] = arr[j][i];
                arr[j][i] = tmp;
            }
        }
    }

    public static void rev(int[] arr) {
//        Two pointer swapping
        int len = arr.length;
        int i = 0, j = len - 1;
        while (i < j) {
            swap(arr, i, j);
            i++;
            j--;
        }

    }

    //    SPIRAL TRAVERSAL OF A MATRIX
    public static void print_spiral(int[][] arr) {
        //T.C --> O(N*M)
        int n = arr.length, m = arr[0].length;
        int tot_ele = n * m;
        int top = 0, right = m - 1, bottom = n - 1, left = 0;
        List<Integer> ans = new ArrayList<>();
        //Go right
        while (top <= bottom && left <= right) {
            //GO RIGHT
            for (int i = left; i <= right; i++) {
                ans.add(arr[top][i]);
                top++;
            }
            //go bottom
            for (int i = top; i <= bottom; i++) {
                ans.add(arr[i][right]);
                right--;
            }

            //go left
            if (top <= bottom) { //edge case what if it is single line
                for (int i = right; i >= left; i--) {
                    ans.add(arr[bottom][i]);
                    bottom--;
                }
            }
            //go top
            for (int i = bottom; i >= top; i--) {
                ans.add(arr[i][left]);
                left++;
            }

        }
    }

    /*Pascal Triangle --> Types of problems
    1. Given R & C print the element      --> ncr
    2. print nth row               --> call function ncr n times for( col 1 to n)  T.C--> N*R
                                    --> pascal_optimal generate row function
    3. given n print whole pascal triangle till n   --> call generate row function n times

    Formula-->   nCr  --> (r-1)C(c-1)
    nCr = n!/(r! * (n-r)!    //too much time complexity
    Shortcut --> n! will always goes till r! as rest gets divided by (n-r)!
     */
    //Brute Force
    public static int ncr(int n, int r) {
        //T.C--> O(R)  S.C --> O(1)
        int res = 1;
        for (int i = 0; i < r; i++) {
            res *= (n - i);
            res = res / (i + 1);
        }

        return res;
    }

    //optimal  --> T.C = O(N)
    public static int pascal_optimal(int n) {
//Observation --> we will just divide the answer by column
        int ans = 1;
        for (int col = 1; col < n; col++) {
            System.out.print(ans + " ");
            ans *= (n - col);
            ans /= col;
        }
        System.out.println(1);
        return ans;
    }

    //    Majority element II
    public List<Integer> majorityElement(int[] arr) {
        //        ele's  > [n/3] times
//Observation --> the answer wil be having at max 2 integers
        ArrayList<Integer> ar = new ArrayList<>();
        int cnt1 = 0, cnt2 = 0, n = arr.length, ele1 = Integer.MIN_VALUE, ele2 = Integer.MIN_VALUE;
        for (int i = 0; i < n; i++) {
            if (cnt1 == 0 && arr[i] != ele2) {
                ele1 = arr[i];
                cnt1++;
            } else if (cnt2 == 0 && arr[i] != ele1) {
                ele2 = arr[i];
                cnt2++;
            } else if (ele1 == arr[i]) {
                cnt1++;
            } else if (ele2 == arr[i]) {
                cnt2++;
            } else {
                cnt1--;
                cnt2--;
            }
        }

        //manual check
        cnt1 = 0;
        cnt2 = 0;
        for (int i = 0; i < n; i++) {
            if (arr[i] == ele1) cnt1++;
            if (arr[i] == ele2) cnt2++;
        }
        if (cnt1 > (n / 3)) {
            ar.add(ele1);
        }
        if (cnt2 > (n / 3)) {
            ar.add(ele2);
        }
        return ar;
    }

    /*  3 SUM
    Brute force --> Run 3 loops and then sort & store in set to check every triplet is unique
    T.C--> O(N^3) * LOG (NO. OF UNIQUE TRIPLETS)
    S.C--> 2 * O(NO. OF UNIQUE TRIPLETS)

    BETTER-> (Hashing) --> Removes 3rd loop

    Optimal-> get rid of this set
     */
    public static List<List<Integer>> Three_Sum_Brute(int[] arr) {
//    return triplets whose sum equals to 0
        List<List<Integer>> ans = new ArrayList<>();
        Set<List<Integer>> st = new HashSet<>();
        int n = arr.length;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                for (int k = j + 1; k < n; k++) {
                    if (arr[i] + arr[j] + arr[k] == 0) {
                        List<Integer> ls = new ArrayList<>();
                        ls.add(arr[i]);
                        ls.add(arr[j]);
                        ls.add(arr[k]);
                        Collections.sort(ls);
                        st.add(ls); //maintains Uniqueness
                    }
                }
            }
        }

        for (List<Integer> ls : st) {
            ans.add(ls);
        }

        return ans;
    }

    public static List<List<Integer>> Three_Sum_Better(int[] arr) {
//    arr[k] = -(arr[i] + arr[j])
//        t.c --> O(N^2) * LOG(M)   , S.C --> O(N) + 2 * O(NO. OF UNIQUE TRIPLETS)

        List<List<Integer>> ans = new ArrayList<>();
        Set<List<Integer>> st = new HashSet<>();

        int n = arr.length;
        for (int i = 0; i < n; i++) {
            Set<Integer> chck = new HashSet<>();
            for (int j = i + 1; j < n; j++) {
                if (chck.contains(-(arr[i] + arr[j]))) {
                    List<Integer> ls = new ArrayList<>();
                    ls.add(arr[i]);
                    ls.add(arr[j]);
                    ls.add(-(arr[i] + arr[j]));
                    Collections.sort(ls);
                    st.add(ls); //maintains Uniqueness
                }
                chck.add(arr[j]);
            }
        }
        for (List<Integer> ls : st) {
            ans.add(ls);
        }
        return ans;
    }

    public static List<List<Integer>> Three_Sum_Optimal(int[] arr) {
//    arr[k] = -(arr[i] + arr[j])
//        t.c --> O(N LOG N)  + (N^2)    , S.C --> O(no. of unique triplets ) for storing triplets
/*Observation --> Eliminate extra space which is used to filter out duplicate triplets
  sort the array instead and formulate algorithm
  TWO POINTER APPROACH 
 */

        List<List<Integer>> ans = new ArrayList<>();
        Arrays.sort(arr);
        int n = arr.length;
        for (int i = 0; i < n; i++) {
            if (i > 0 && arr[i] == arr[i - 1]) continue;
            int j = i + 1, k = n - 1;
            while (j < k) {
                int sum = arr[i] + arr[j] + arr[k];
                if (sum < 0) {
                    j++;
                } else if (sum > 0) {
                    k--;
                } else {
                    List<Integer> tmp = new ArrayList<>();
                    tmp.add(arr[i]);
                    tmp.add(arr[j]);
                    tmp.add(arr[k]);
                    ans.add(tmp);
                    j++;
                    k--;
                    while (j < k && arr[j] == arr[j - 1]) j++;
                    while (j < k && arr[k] == arr[k + 1]) k--;
                }

            }
        }

        return ans;
    }

    public static List<List<Integer>> Four_Sum_Optimal(int[] arr, int target) {
        //T.C --> O(NLOG N) + O(N^3)
        List<List<Integer>> ans = new ArrayList<>();
        Arrays.sort(arr);
        int n = arr.length;
        for (int i = 0; i < n; i++) {
            if (i > 0 && arr[i] == arr[i - 1]) continue;
            for (int j = i + 1; j < n; j++) {
                if (j > i + 1 && arr[j] == arr[j - 1]) continue;
                //two pointer's
                int k = j + 1, l = n - 1;
                while (k < l) {
                    long sum = (long) arr[i] + arr[j] + arr[k] + arr[l];
                    if (sum < target) {
                        k++;
                    } else if (sum > target) {
                        l--;
                    } else {
                        List<Integer> tmp = new ArrayList<>();
                        tmp.add(arr[i]);
                        tmp.add(arr[j]);
                        tmp.add(arr[k]);
                        tmp.add(arr[l]);
                        ans.add(tmp);
                        k++;
                        l--;
                        //uniqueness
                        while (k < l && arr[k] == arr[k - 1]) k++;
                        while (k < l && arr[l] == arr[l + 1]) l--;
                    }
                }
            }
        }
        return ans;
    }

    //   Count   Subarray with xor k
    //optimal
    public static int cnt_subarr_xor_k(int[] arr, int k) {
        int cnt = 0;
        HashMap<Integer, Integer> map = new HashMap<>();
        map.put(0, 1);
        int prexor = 0;
        for (int ele : arr) {
            prexor ^= ele;
            if (map.containsKey(prexor ^ k)) {
                cnt += map.get(prexor ^ k);
            }
            map.put(prexor, map.getOrDefault(prexor, 0) + 1);
        }
        return cnt;
    }

    //    Merge Intervals --optimized
    public static int[][] merge_Intervals(int[][] arr) {
//        [[1,4],[2,3]]
//        T.C--> O(NLOGN) + O(N)
//        S.C--> O(N)
        sort_intervals(arr);
        List<int[]> ls = new ArrayList<>();
        int n = arr.length;

        for (int i = 0; i < n; i++) {
            int lst = ls.get(ls.size() - 1)[1];
            if (ls.isEmpty() || lst < arr[i][0]) {
                //add interval
                ls.add(arr[i]);
            } else {
                ls.get(ls.size() - 1)[1] = Math.max(arr[i][1], ls.get(ls.size() - 1)[1]);
            }
        }

        return ls.toArray(new int[ls.size()][]);
    }

    public static void sort_intervals(int[][] arr) {
        Arrays.sort(arr, (a, b) -> {
            if (a[0] == b[0]) {
                return Integer.compare(a[1], b[1]);
            } else {
                return Integer.compare(a[0], b[0]);
            }
        });
    }

    //MERGE 2 SORTED ARRAYS WITHOUT USING EXTRA SPACE
//Brute Force

    public static void merge_arr_ws(int[] arr1, int[] arr2) {
        //T.C --> O(N+M)  + O(N+M)
//        S.C --> O(N+M)
        int n = arr1.length, m = arr2.length;
        int[] arr = new int[n + m];
        int left = 0, right = 0, index = 0;
        while (left < n && right < m) {
            if (arr1[left] <= arr2[right]) {
                arr[index++] = arr1[left++];
            } else {
                arr[index++] = arr2[right++];
            }
        }
        while (left < n) {
            arr[index++] = arr1[left++];
        }
        while (right < m) {
            arr[index++] = arr2[right++];
        }
//fill arrays
        for (int i = 0; i < n + m; i++) {
            if (i < n) {
                arr1[i] = arr[i];
            } else {
                arr2[i - n] = arr[i];
            }
        }
    }

    //optimal 1--> Two pointer
    public static void merge_arr_ws_optimal1(int[] arr1, int[] arr2) {
        //T.C --> O(MIN(N,M) + O(NLOGN) + O(MLOGM)
        int n = arr1.length, m = arr2.length;
        int left = n - 1, right = 0;
        //swap until breakpoint
        while (left >= 0 && right < m) {
            if (arr1[left] >= arr2[right]) {
                swap1(arr1, arr2, left, right);
                left--;
                right++;
            } else {
                //everything is at its correct place from here
                break;
            }
        }
        Arrays.sort(arr1);
        Arrays.sort(arr2);
    }

    public static void swap1(int[] arr1, int[] arr2, int left, int right) {
        if (arr1[left] > arr2[right]) {
            int tmp = arr1[left];
            arr1[left] = arr2[right];
            arr2[right] = tmp;
        }
    }

    //Gap Algorithm
    public static void merge_arr_ws_optimal2(int[] arr1, int[] arr2) {
//  gap  =  (n+m)/2
        //t.c--> Olog2(n+m) + o(n+m)
        int n = arr1.length, m = arr2.length, len = n + m;
        int gap = (len) / 2 + len % 2;  //ceil
        while (gap > 0) {
            int left = 0;
            int right = gap;
            while (right < len) {
                //pointer lies in both
                if (left < n && right >= n) {
                    swap1(arr1, arr2, left, right - n);
                } else if (left < n && right < n) { //arr1
                    swap1(arr1, arr1, left, right);
                } else {  //arr2
                    swap1(arr2, arr2, left - n, right - n);
                }
                left++;
                right++;
            }
            if (gap == 1) break;
            gap = (gap / 2) + (gap % 2);
        }
    }
//find the repeating and missing number
//two approaches maths and xor


    //Approach 1--> Maths
    public static int[] find_rm(int[] arr) {
        //T.C-- > O(N)
        int r = -1, m = -1, n = arr.length, s = 0, sq = 0;
        for (int i = 0; i < n; i++) {
            s += arr[i];
            sq += (arr[i] * arr[i]);
        }
//equations
        int sn = n * (n + 1) / 2;
        int sn2 = n * (n + 1) * (2 * n + 1) / 6;
        //  x-y = ele
        int xMy = s - sn;
        int xPy = (sq - sn2) / xMy;
        r = (xMy + xPy) / 2;
        m = r - xMy;

        return new int[]{r, m};
    }

    //Approach 2 --> xor
    public static int[] find_rm_xor(int[] arr) {
        int n = arr.length, xr = 0;
        //xor array and n natural number
        for (int i = 0; i < n; i++) {
            xr ^= arr[i];
            xr ^= (i + 1);
        }
//        Intuition --> since they are 2 different numbers for sure means
//        they are bound to be different at one bit
//find rightmost set bit number -position
//        int bitNo = 0;
//        while (true) {
//            //shift bit until its set
//            if ((xr & (1 << bitNo)) != 0) {
//                break;
//            }
//            bitNo++;
//        }
        int number = xr & -xr;  //not position its a whole value of that bit
        // find repeating and missing
        int one = 0;
        int zero = 0;

        for (int i = 0; i < n; i++) {
            //part of 1 club
//            if ((arr[i] & (1 << bitNo)) != 0) {
            if ((arr[i] & number) != 0) {
                one ^= arr[i];
            } else {
                zero ^= arr[i];
            }
            //also do for natural numbers
//            if (((i + 1) & (1 << bitNo)) != 0) {
            if (((i + 1) & number) != 0) {
                one ^= (i + 1);
            } else {
                zero ^= (i + 1);
            }
        }

        int cnt = 0;
        for (int i = 0; i < n; i++) {
            if (one == arr[i]) cnt++;
        }
        if (cnt == 2) return new int[]{one, zero};

        return new int[]{zero, one};
    }

    //Count Inversions  -->  i<j , arr[i]> arr[j]

    //    Intuition --> Merge Sort  (split --> sort --> merge )
//    integrate this logic (Given two sorted array cnt inversions) in merge sort
    //merge two sorted arrays
    public static int merge(int[] arr, int low, int mid, int high) {
//        t.c --> O(NLOGN) , S.C --> O(n)
        // Merges two sorted arrays and counts inversions
        int cnt = 0;
        int i = low, j = mid + 1, k = 0;
        int[] res = new int[high - low + 1];  // Temporary array of the required size

        //process of merging
        while (i <= mid && j <= high) {
            if (arr[i] <= arr[j]) {
                res[k] = arr[i];
                i++;
            } else {
                res[k] = arr[j];
//                cnt += (mid - i + 1);  // All elements from i to mid are greater than arr[j]
                j++;
            }
            k++;
        }

        // Copy remaining elements of left subarray, if any
        while (i <= mid) {
            res[k] = arr[i];
            i++;
            k++;
        }

        // Copy remaining elements of right subarray, if any
        while (j <= high) {
            res[k] = arr[j];
            j++;
            k++;
        }

        // Copy the merged elements back into the original array
        for (int l = 0; l < res.length; l++) {
            arr[low + l] = res[l];
        }

        return cnt;
    }

    // Recursive merge sort and inversion count function
    public static int merge_sort(int[] arr, int left, int right) {
        if (left >= right) return 0;  // Base case: If the array has one or no elements, no inversions

        int cnt = 0;
        int mid = left + (right - left) / 2;

        // Count inversions in the left half
        cnt += merge_sort(arr, left, mid);

        // Count inversions in the right half
        cnt += merge_sort(arr, mid + 1, right);
        // Count inversions during the merge step
        cnt += merge(arr, left, mid, right);

        return cnt;
    }


    // Function to count inversions in the array
    public static int cnt_inversion(int[] arr) {
        return merge_sort(arr, 0, arr.length - 1);
    }

    //    Reverse Pairs
//    pair is only possible when something is from the left array and something is from the right
    public int reversePairs(int[] arr) {
        return merge_sort_rev_pairs(arr, 0, arr.length - 1);
    }

    public static int count_pairs(int arr[], int low, int mid, int high) {
        //process to count reverse pairs
// Approach --> slight change in two pointer technique used in cnt inversion
        //implementation
        int cnt = 0;
        int y = mid + 1;
        for (int x = low; x <= mid; x++) {
            while (y <= high && (long) arr[x] > (long) 2 * arr[y]) {
                y++;
            }
            cnt += (y - (mid + 1));
        }
        return cnt;
    }

    // Recursive merge sort and inversion count function
    public static int merge_sort_rev_pairs(int[] arr, int left, int right) {
        if (left >= right) return 0;  // Base case: If the array has one or no elements, no inversions

        int cnt = 0;
        int mid = left + (right - left) / 2;

        // Count inversions in the left half
        cnt += merge_sort(arr, left, mid);

        // Count inversions in the right half
        cnt += merge_sort(arr, mid + 1, right);
        cnt += count_pairs(arr, left, mid, right);
        // Count inversions during the merge step
        merge(arr, left, mid, right);

        return cnt;
    }

    //Maximum product subarray
    public static int max_prod_subArr(int[] arr) {
        /*OBSERVATION --> Either the answer will be somewhere at the prefix or somewhere at the suffix
        --> array has all +ve
        --> array has even -ve's
        --> array has odd +-ve's
        --> array has zeroes
         */
        //T.C --> O(N)
        int ans = Integer.MIN_VALUE, n = arr.length;
        int prefx = 1, suffx = 1;
        for (int i = 0; i < n; i++) {
            if (prefx == 0) prefx = 1;
            if (suffx == 0) suffx = 1;
            prefx *= arr[i];
            suffx *= arr[n - i - 1];
            ans = Math.max(ans, Math.max(prefx, suffx));
        }
        return ans;
    }

    /*BINARY SEARCH PATTERN
    1.Typical bs on arr
    2.Binary search on answers min or max
    3. bs on min(max) or max(min)


     */
    // Iterative  Binary Search
    public static boolean bin_sear(int[] arr, int ele) {
        int n = arr.length;
        int i = 0, j = n - 1;
        while (i <= j) {
            int mid = i + (j - i) / 2;
            if (arr[mid] == ele) return true;
            else if (arr[mid] < ele) {
                i = mid + 1;
            } else {
                j = mid - 1;
            }
        }
        return false;
    }

    public static boolean bin_sear_rec(int[] arr, int ele, int i, int j) {
//base case
        if (i > j) return false;
        int mid = i + (j - i) / 2;
        if (arr[mid] == ele) return true;
        if (arr[mid] < ele) {
            return bin_sear_rec(arr, ele, mid + 1, j);
        }
        if (arr[mid] > ele) {
            return bin_sear_rec(arr, ele, i, mid - 1);
        }

        return false;
    }

    //    Lower Bound and Upper bound
    public static int lower_bound(int[] arr, int ele) {
//        arr[idx] >=ele  //the smallest idx
//        if there is noone then includes last hypothetical idx that is size of the array
        int n = arr.length;
        int idx = n;
        int i = 0, j = n - 1;
        while (i <= j) {
            int mid = i + (j - i) / 2;
            if (arr[mid] >= ele) {
                idx = mid;
                //search for smallest idx
                j = mid - 1;
            } else if (arr[mid] < ele) {
                i = mid + 1;
            }
        }
        return idx;
    }

    public static int Upper_bound(int[] arr, int ele) {
//   Strictly Greater
        int n = arr.length;
        int idx = n;
        int i = 0, j = n - 1;
        while (i <= j) {
            int mid = i + (j - i) / 2;
            if (arr[mid] > ele) {
                idx = mid;
                //search for smallest idx
                j = mid - 1;
            } else if (arr[mid] <= ele) {
                i = mid + 1;
            }
        }
        return idx;
    }

    //Floor and ceil in sorted array
    public static int floor_ceil(int[] arr, int ele) {
/*
Floor --> Largest number in arr arr[i]<=x
Ceil --> smallest number in arr arr[i]>=x   --> lower bound
 */

        //for floor
        int n = arr.length;
        int ans = -1;
        int i = 0, j = n - 1;
        while (i <= j) {
            int mid = i + (j - i) / 2;
            if (arr[mid] <= ele) {
                ans = arr[mid];
                //search for largest number
                i = mid + 1;  //go right
            } else if (arr[mid] > ele) {
                j = mid - 1; //go left
            }
        }
        return ans;
    }
    //First & Last occurence also count of an element in an sorted array with duplicates

    public static void occurence(int[] arr, int ele) {
        int first_occ = -1, last_occ = -1;
        int n = arr.length;
        int i = 0, j = n - 1;
        while (i <= j) {
            int mid = i + (j - i) / 2;
            if (arr[mid] == ele) {
                first_occ = mid;
                j = mid - 1;
            } else if (arr[mid] < ele) {
                i = mid + 1;
            } else {
                j = mid - 1;
            }
        }


        i = 0;
        j = n - 1;
        while (i <= j) {
            int mid = i + (j - i) / 2;
            if (arr[mid] == ele) {
                last_occ = mid;
                i = mid + 1;
            } else if (arr[mid] < ele) {
                i = mid + 1;
            } else {
                j = mid - 1;
            }
        }
        int cnt = last_occ - first_occ + 1;
    }

    //Search in a rotated sorted array unique elements only
    /*
    Identify the sorted half --  left/right
     */
    public static int search(int[] arr, int tar) {
        int n = arr.length, i = 0, j = n - 1;

        while (i <= j) {
            int mid = i + (j - i) / 2;
            if (arr[mid] == tar) {
                return mid;
            }
            //identify the sorted
            if (arr[i] <= arr[mid]) {//left half is sorted
                //check if tar lies
                if (tar >= arr[i] && tar < arr[mid]) { //lies
                    j = mid - 1;
                } else { //not lies
                    i = mid + 1;
                }
            } else { //right half is sorted
                //check if tar lies
                if (tar <= arr[j] && tar > arr[mid]) { //lies
                    i = mid + 1;
                } else { //not lies
                    j = mid - 1;
                }
            }
        }
        return -1;
    }

    //    Contains Duplicates
//    the only condition we need to handle is  arr[i]==arr[mid]==arr[j]
    public static boolean searchII(int[] arr, int tar) {
        //worst case T.C --> O(N/2)   array might shrink half of the time if arr contains lots of duplicates

        int n = arr.length, i = 0, j = n - 1;

        while (i <= j) {
            int mid = i + (j - i) / 2;
            if (arr[mid] == tar) {
                return true;
            }
            //edge case
            if (arr[i] == arr[mid] && arr[j] == arr[mid]) {
                i++;
                j--;    //just shrink the array
                continue;
            }
            //identify the sorted
            if (arr[i] <= arr[mid]) {//left half is sorted
                //check if tar lies
                if (tar >= arr[i] && tar < arr[mid]) { //lies
                    j = mid - 1;
                } else { //not lies
                    i = mid + 1;
                }
            } else { //right half is sorted
                //check if tar lies
                if (tar <= arr[j] && tar > arr[mid]) { //lies
                    i = mid + 1;
                } else { //not lies
                    j = mid - 1;
                }
            }
        }
        return false;
    }

    //Find miniumum in rotated and sorted array
    public static int findMin(int[] arr) {
        /*Intuition
         -->identify rotating point as it always be minimum
         --> point will always remains in unsorted part
         --> Discard sorted part but sorted part may or may not have min element
         --> pick min of sorted park before discarding
         */
        int ans = Integer.MAX_VALUE, n = arr.length, i = 0, j = n - 1;
        while (i <= j) {
            int mid = i + (j - i) / 2;
            //if the search space is already sorted
            if (arr[i] <= arr[j]) return arr[i];
            if (arr[mid] >= arr[i]) {
                //sorted half
                ans = Math.min(ans, arr[i]);
                //discard
                i = mid + 1;
            } else {
                ans = Math.min(ans, arr[mid]);
                j = mid - 1;
            }
        }
        return ans;
    }

    //find out how many times array has been rotated
    public static int cnt_rot(int[] arr) {
        int ans = Integer.MAX_VALUE;
        int idx = 0, n = arr.length, i = 0, j = n - 1;
        while (i <= j) {
            int mid = i + (j - i) / 2;
            //if the search space is already sorted
            if (arr[i] <= arr[j]) {
                return i;
            }
            if (arr[mid] >= arr[i]) {
                //sorted half
                if (arr[i] < ans) {
                    ans = arr[i];
                    idx = i;
                }
                //discard
                i = mid + 1;
            } else {
                if (arr[mid] < ans) {
                    ans = arr[mid];
                    idx = mid;
                }
                j = mid - 1;
            }
        }
        return idx;

        //for duplicates just shrink the array
    }
    //Single element in a sorted array

    public static int singleNonDuplicate(int[] arr) {
        /*
        Observation --> even odd approach
        (even,odd) --> ele lies on right half
        (odd,even) --> ele lies on left half

         */
        int n = arr.length, i = 1, j = n - 2;
        //edge cases
        if (n == 1 || arr[0] != arr[1]) return arr[0];
        if (arr[n - 1] != arr[n - 2]) return arr[n - 1];
        while (i <= j) {
            int mid = i + (j - i) / 2;
            if (arr[mid] != arr[mid + 1] && arr[mid] != arr[mid - 1]) return arr[mid];
            //standing at even idx
            if (mid % 2 == 0 && arr[mid] == arr[mid + 1] || mid % 2 == 1 && arr[mid] == arr[mid - 1]) {
                //ele lies on right half
                i = mid + 1;
            } else {
                j = mid - 1;
            }
        }
        return -1;
    }

    //Find Peak Element
    //brute force
    public int findPeakElement(int[] arr) {
        int n = arr.length;
        for (int i = 0; i < n; i++) {
            if ((i == 0 || arr[i] > arr[i - 1]) && (i == n - 1 || arr[i] > arr[i + 1])) {
                return i;
            }
        }
        return 0;
    }

    public int findPeakElement_optimal(int[] arr) {
//peak will always remain above increasing curve  --> right side
// mid can be peak
//On decreasing curve peak will remain in --->  left side
        int n = arr.length, i = 0, j = n - 1;
        if (n == 1) return 0;
        while (i <= j) {
            int mid = i + (j - i) / 2;
            //mid is peak
            if ((mid == 0 || arr[mid] > arr[mid - 1]) && (mid == n - 1 || arr[mid] > arr[mid + 1])) {
                return mid;
            }
            //mid lies on increasing half
            if ((mid == 0 || arr[mid] > arr[mid - 1]) && (mid == n - 1 || arr[mid] < arr[mid + 1])) {
                i = mid + 1;
            }
            //mid lies on decreasing half
            else {
                j = mid - 1;
            }
        }

        return -1;
    }

    //Find the sqrt of an integer using Binary Search
    public static int sqrt(int n) {
        int i = 0, j = n, ans = 0;
        while (i <= j) {
            int mid = i + (j - i) / 2;
            if (mid * mid < n) { //floor
                ans = mid;
                i = mid + 1;
            } else {
                j = mid - 1;
            }
        }

        return ans;
    }

    //Find the Nth root of an Integer
    public static int Nth_root(int x, int n) {
        //x = 27   n = 3   ==> 3 * 3 * 3 = 27
        //Brute force  T.C--- O(X)*O(LOG2N)
        for (int i = 0; i <= x; i++) {
            int tot = (int) myPow(i, n);
            if (tot == x) {
                return i;
            } else if (tot > n) {
                return -1;
            }

        }
        return -1;
    }

    //Using Binary Search
    public static int Nth_root_optimal(int x, int n) {
        int i = 1, j = x;
        while (i <= j) {
            int mid = i + (j - i) / 2;
//            int tot = (int) myPow(mid, n);
            int MidN = MidN(mid, x, n);
            if (MidN == 1) return mid;
            else if (MidN == 0) {
                i = mid + 1;
            } else {
                j = mid - 1;
            }

        }

        return -1;
    }
    //Optimized myPow functio for this problem only to prevent overflow

    public static int MidN(int mid, int x, int n) {
        long ans = 1;
        for (int i = 1; i < n; i++) {
            ans *= ans;
            if (ans > x) return 2;
        }
        if (ans == x) return 1;

        return 0; //ans is < than x
    }

    public static double myPow(double x, int n) {
        /*Also Handle if n is -ve
        if(n<0){
            x = 1/x;
            n = -n;
            return myPow1(x,n);
        }
         */

        // BY LOOP  --> O(N)
        //Using power exponential method  = o(log2n)
 /*
 2^21  -->  odd --> 2 * 2^20
 2^20  -->  even  --> 2 * (2^2)20/2  means 4^10
We just dropped of by half that's the intuition
  */
        if (n == 0) return 1;
        if (n % 2 == 0) { //even
            return myPow(x * x, n / 2);
        } else { //odd
            return x * myPow(x * x, n / 2);
        }
    }

    //KOKO Eating bananas

    //Brute Force
    public int minEatingSpeed(int[] piles, int h) {
        int max = Integer.MIN_VALUE;
        for (int ele : piles) {
            max = Math.max(ele, max);
        }
        for (int bPh = 1; bPh <= max; bPh++) {
            //figure out the total time
            long tot_time = 0;
            for (int i = 0; i < piles.length; i++) {
                tot_time += (int) Math.ceil((double) piles[i] / bPh); //ceil

            }
            if (tot_time <= h) return bPh;
        }
        return -1;
    }

    //Koko Optimized Using Binary Search
    public int minEatingSpeed_optimized(int[] piles, int h) {
        int max = Integer.MIN_VALUE;
        for (int ele : piles) {
            max = Math.max(ele, max); // Find the maximum pile
        }

        int ans = max;
        int i = 1, j = max;
        while (i <= j) {
            int mid = i + (j - i) / 2;
            if (canFinishInTime(piles, mid, h)) {
                ans = mid; // Store the possible answer (minimum speed)
                j = mid - 1; // Try for a smaller speed
            } else {
                i = mid + 1; // Increase speed
            }
        }
        return ans; // Return the minimum eating speed
    }

    private boolean canFinishInTime(int[] piles, int speed, int h) {
        long totalHours = 0; // Use long to avoid overflow
        for (int ele : piles) {
            totalHours += (int) Math.ceil((double) ele / speed); // Calculate total time
        }
        return totalHours <= h; // Return true if Koko can finish in h hours
    }

    //Minimum Number of days to make M bouquets
    public int minDays(int[] bloomDay, int bouq, int flow) {
        //T.C --> O(N) *O(LOG2N)
        int n = bloomDay.length, max = bloomDay[0];
        if (n < flow * bouq) return -1;
        int min = bloomDay[0];
        for (int ele : bloomDay) {
            max = Math.max(ele, max);
            min = Math.min(ele, min);

        }
        int i = min;
        int minDays = max, j = max;

        while (i <= j) {
            int days = i + (j - i) / 2;
            if (canBloom(bloomDay, bouq, days, flow)) {
                minDays = days;
                j = days - 1;
            } else {
                i = days + 1;
            }

        }
        return minDays;
    }

    private boolean canBloom(int[] bloomDay, int bouq, int days, int flow) {

        boolean adjacent = true;
        int cur = 0;
        for (int i = 0; i < bloomDay.length; i++) {
            if (bloomDay[i] <= days) {
                cur++;
            } else {
                //start new
                cur = 0;
            }
            if (cur == flow) {
                bouq--;
                cur = 0;
            }
            if (bouq == 0) {
                return true;
            }
        }
        return false;
    }

    //    Find the Smallest Divisor Given a Threshold
    public int smallestDivisor(int[] arr, int threshold) {
        int i = 1, j = threshold, ans = threshold;
        while (i <= j) {
            int div = i + (j - i) / 2;
            if (possible(arr, div, threshold)) {
                ans = div;
                j = div - 1;
            } else {
                i = div + 1;
            }
        }
        return ans;
    }

    private boolean possible(int[] arr, int div, int threshold) {
        int res = 0;
        for (int ele : arr) {
            res += (int) Math.ceil((double) ele / (double) div);
        }
        return res <= threshold;
    }

    //Capacity to ship packages within D days
    public int shipWithinDays(int[] weights, int days) {
        int j = weights[0], sum = 0;
        for (int weight : weights) {
            j = Math.max(j, weight);
            sum += weight;
        }
        int i = j;
        int cap = i;
        j = sum;
        while (i <= j) {
            int leastCap = i + (j - i) / 2;
            if (possible_capacity(weights, days, leastCap)) {
                cap = leastCap;
                j = leastCap - 1;
            } else {
                i = leastCap + 1;
            }

        }
        return cap;
    }

    private boolean possible_capacity(int[] weights, int days, int leastCap) {
        int cur_cap = leastCap;
        //within
        int daysNeeded = 1;
        for (int weight : weights) {

            if (cur_cap >= weight) {
                cur_cap -= weight;
            } else {
                daysNeeded++;
                cur_cap = leastCap - weight;
            }
            if (daysNeeded > days) return false;
        }
        return true;
    }
    //kth Missing Number

    //Brute Force
    public int findKth_Positive(int[] arr, int k) {
        for (int ele : arr) {
            if (ele > k) return k;
            else k++;
        }
        return k;
    }

    //Optimized using Binary Search
    public int findKthPositive(int[] arr, int k) {
        //Approach --> get two nearby indexes
        //Discard left /right on the basis of  missing numbers and k

        int n = arr.length, i = 0, j = n - 1;
        while (i <= j) {
            int mid = i + (j - i) / 2;
            int missing_num = arr[mid] - (mid + 1);
            if (missing_num < k) {
                i = mid + 1;
            } else {
                j = mid - 1;
            }
        }
//Opposite polarity --> At some point high will cross low and we get range of two nearby idx's
        //derivation

//        arr[j] + (k - (arr[j] - (j+1))) ==   j+1 +k
        //j +1 is also equals to i (as i is just ahead)
        return i + k;
    }

    //    Aggressive cows
//   (min distance b/w cows is maximum)
    public static int Aggressive_cows(int[] arr, int cows) {
        //range
        Arrays.sort(arr); //sort coordinates
        int n = arr.length, i = 1, j = arr[n - 1] - arr[0];
        while (i <= j) {
            int mid = i + (j - i) / 2;
            if (CanWePlace(arr, mid, cows)) {
                i = mid + 1;
            } else {
                j = mid - 1;
            }

        }
        //opposite polarity means not possible pointer will become the answer always

        return j;  //as per polarity high will always remain the answer
    }

    private static boolean CanWePlace(int[] arr, int dist, int cows) {

        int tot_cows_placed = 1, lastCow = arr[0];
        for (int i = 1; i < arr.length; i++) {
            if (arr[i] - lastCow >= dist) {
                tot_cows_placed++;
                lastCow = arr[i];
            }
            if (tot_cows_placed == cows) return true;
        }
        return false;
    }

    //Book Allocation or Painters partition leetcode or split array largest sum min(max)
    //maximum number of pages allocated is minimum
    public static int Book_Allocation(int[] pages, int stu) {
        int n = pages.length;
        if (n < stu) return -1; //each student should have atleast 1 book

//try later
        return -1;
    }

    //Minimise max distance between gas stations
    //Brute Force
    public static int min_max_dist_btwn_gas_station(int[] arr, int k) {
        //T.C --> K * N  + N
        // {1, 13, 17, 23}
        int n = arr.length;
        int[] howMany = new int[n - 1];  //keeps track of how many have been placed
        for (int gasStation = 1; gasStation <= k; gasStation++) {
//place at the max value sector so keeps the track of it
            int maxLen = -1, maxIdx = -1;
            for (int i = 0; i < n - 1; i++) {
                //get max difference sector to place
                int diff = arr[i + 1] - arr[i];
                int section_length = diff / (howMany[i] + 1);
                if (section_length > maxLen) {
                    maxLen = section_length;
                    maxIdx = i;
                }
            }
            howMany[maxIdx]++;
        }

        int maxAns = -1;
        for (int i = 0; i < n - 1; i++) {
            int diff = arr[i + 1] - arr[i];
            int section_length = diff / (howMany[i] + 1);
            maxAns = Math.max(maxAns, section_length);
        }
        return maxAns;
    }

    //OPTIMIZED
    public static int min_max_dist_btwn_gas_station_optimized(int[] arr, int k) {
        //T.C -->  NLOGN  + KLOGN  , S.C -->   O(N-1)
        // {1, 13, 17, 23}
        int n = arr.length;
        PriorityQueue<Pair> pq = new PriorityQueue<>((a, b) -> b.first - a.first); //max-heap

        int[] howMany = new int[n - 1];  //keeps track of how many have been placed

        for (int i = 0; i < n - 1; i++) {
            //get max difference sector to place
            int diff = arr[i + 1] - arr[i];
            pq.add(new Pair(diff, i));
        }
        for (int gasStation = 1; gasStation <= k; gasStation++) {
//place at the max value sector so keeps the track of it
            Pair tp = pq.poll();
            int secIdx = tp.second;
            howMany[secIdx]++;
            //update
            int diff = arr[secIdx + 1] - arr[secIdx];
            int section_length = diff / (howMany[secIdx] + 1);
            pq.add(new Pair(section_length, secIdx));

        }
        return pq.peek().first;
    }

    static class Pair {
        int first;
        int second;

        public Pair(int first, int second) {
            this.first = first;
            this.second = second;
        }
    }

    //Median of two sorted arrays
    /*Brute Force
   merge them and if n is even then   [n]/2 +[n/2-1]
   if odd simple   [n/2] is ans
   */
    //Optimizing space  using hypothetical tracker
    public static double median_arr(int[] arr1, int[] arr2) {
        int n = arr1.length, m = arr2.length, t = n + m;
        //needed indexes
        int idx2 = t / 2, idx1 = idx2 - 1;  //required for even
        int ele1 = -1, ele2 = -1;
        //for odd only required   idx2 and ele2
        int cnt = 0, i = 0, j = 0;
        while (i < n && j < m) {
            if (arr1[i] <= arr2[j]) {
                if (cnt == idx1) ele1 = arr1[i];
                if (cnt == idx2) ele2 = arr1[i];
                cnt++;
                i++;
            } else {
                if (cnt == idx1) ele1 = arr2[j];
                if (cnt == idx2) ele2 = arr2[j];
                cnt++;
                j++;
            }
        }
        while (i < n) {
            if (cnt == idx1) ele1 = arr1[i];
            if (cnt == idx2) ele2 = arr1[i];
            cnt++;
            i++;
        }
        while (j < m) {
            if (cnt == idx1) ele1 = arr2[j];
            if (cnt == idx2) ele2 = arr2[j];
            cnt++;
            j++;
        }
        if (t % 2 == 1) return ele2;

        return ((double) (ele1 + ele2) / 2.0);
    }
    //Approach 3 using binary search and split technique

    //TRY KTH ELEMENT OF TWO SORTED ARRAYS

    //    Find  the row with maximum 1's
//    public static int[] Max_ones(int[][] arr) {
//      T.C --> On LOG M
//
//    }
    /*Search element in a 2D matrix
    flatten 2d array hypothetically and apply binary search
    convert the 1D convert ---> 2D coordinate

    row  = idx/m
    col = idx%m
T.C --> OLOG2(N*M)
     */

}

